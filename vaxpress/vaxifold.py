from aio_pika import connect, Message
import json
import uuid
import asyncio
from struct import unpack
import numpy as np
from scipy import sparse  # type: ignore
from tqdm import tqdm  # type: ignore
from tqdm.asyncio import tqdm_asyncio as atqdm  # type: ignore

from concurrent.futures import Executor, ProcessPoolExecutor
from importlib.metadata import version

from .log import log


class LocalVaxiFold:

    def __init__(self, folding_engine: str, partition_engine: str):
        self.folding_engine = folding_engine
        self.partition_engine = partition_engine

    async def __aenter__(self):
        return self

    async def __aexit__(self, exc_type, exc_value, traceback):
        pass

    async def fold(self, seqs, executor: Executor = None):  # type: ignore[assignment]
        executor = executor or ProcessPoolExecutor()
        if self.folding_engine == "viennarna":
            func = self.viennarna_fold
        elif self.folding_engine == "linearfold":
            func = self.linearfold
        else:
            raise ValueError(f"Unknown engine: {self.folding_engine}")

        total_tasks = len(seqs)
        progress_bar = atqdm(total=total_tasks, desc="Folding seqs", unit="requests")

        async def fold_with_progress(seq):
            result = await asyncio.wrap_future(executor.submit(func, seq))
            progress_bar.update(1)
            return result

        tasks = [fold_with_progress(seq) for seq in seqs]
        results = await asyncio.gather(*tasks)

        progress_bar.close()
        return results

    async def partition(self, seqs, executor: Executor = None):  # type: ignore[assignment]
        executor = executor or ProcessPoolExecutor()
        if self.partition_engine == "viennarna":
            func = self.viennarna_partition
        elif self.partition_engine == "linearpartition":
            func = self.linearpartition
        else:
            raise ValueError(f"Unknown engine: {self.partition_engine}")

        total_tasks = len(seqs)
        progress_bar = atqdm(total=total_tasks, desc="Partition seqs", unit="requests")

        async def partition_with_progress(seq):
            result = await asyncio.wrap_future(executor.submit(func, seq))
            progress_bar.update(1)
            return result

        tasks = [partition_with_progress(seq) for seq in seqs]
        results = await asyncio.gather(*tasks)

        progress_bar.close()
        return results

    @staticmethod
    def viennarna_fold(seq):
        try:
            import ViennaRNA as RNA  # type: ignore
        except ImportError:
            raise ImportError(
                "ViennaRNA is not installed. \n try: `pip install ViennaRNA`"
            )
        assert version("ViennaRNA") >= "2.4.0", "ViennaRNA version >= 2.4.0 is required"

        struct, mfe = RNA.fold(seq)
        return {"structure": struct, "free_energy": mfe}

    @staticmethod
    def linearfold(seq):
        try:
            import linearfold  # type: ignore
        except ImportError:
            raise ImportError(
                "LinearFold is not installed. \n try: `pip install linearfold-unofficial`"
            )
        assert (
            version("linearfold-unofficial") >= "0.1"
        ), "LinearFold version >= 0.1 is required"

        struct, mfe = linearfold.fold(seq)
        return {"structure": struct, "free_energy": mfe}

    @staticmethod
    def viennarna_partition(seq):
        try:
            import ViennaRNA as RNA  # type: ignore
        except ImportError:
            raise ImportError(
                "ViennaRNA is not installed. \n try: `pip install ViennaRNA`"
            )
        assert version("ViennaRNA") >= "2.4.0", "ViennaRNA version >= 2.4.0 is required"

        fc = RNA.fold_compound(seq)
        fold, fe = fc.pf()
        bpp = np.array(fc.bpp())[1:, 1:]
        bp_dtype = [("i", "i4"), ("j", "i4"), ("prob", "f8")]
        indices = np.nonzero(bpp)
        probabilities = bpp[indices]
        bpmtx = np.zeros(probabilities.size, dtype=bp_dtype)
        bpmtx["i"] = indices[0]
        bpmtx["j"] = indices[1]
        bpmtx["prob"] = probabilities
        results = {
            "structure": fold,
            "free_energy": fe,
            "bpp": bpmtx,
            "pi_array": sparse.coo_matrix(bpp + bpp.T),
        }
        assert len(fold) == len(seq), "Partitioning failed for sequence: %s" % seq
        return results

    @staticmethod
    def linearpartition(seq):
        try:
            import linearpartition  # type: ignore
        except ImportError:
            raise ImportError(
                "LinearPartition is not installed. \n try: `pip install linearpartition-unofficial`"
            )
        assert (
            version("linearpartition-unofficial") >= "0.3"
        ), "LinearPartition version >= 0.3 is required"

        pred = linearpartition.partition(seq)
        results = {
            "structure": pred["structure"],
            "free_energy": pred["free_energy"],
            "bpp": pred["bpp"],
        }

        base_pair_i = results["bpp"]["i"]
        base_pair_j = results["bpp"]["j"]
        base_pair_prob = results["bpp"]["prob"]
        pair_prob_one_sided = sparse.coo_matrix(
            (base_pair_prob, (base_pair_i, base_pair_j)), shape=(len(seq), len(seq))
        )
        pair_prob = pair_prob_one_sided + pair_prob_one_sided.T
        results["pi_array"] = pair_prob

        return results


class AsyncVaxiFoldClient:

    def __init__(
        self, host, port, user, passwd, queue, folding_engine, partition_engine
    ):
        self.host = host
        self.port = port
        self.user = user
        self.passwd = passwd
        self.queue = queue
        self.folding_engine = folding_engine
        self.partition_engine = partition_engine
        self.progress_bars = {}
        self.futures = {}
        self.results = {}

    async def connect(self):
        self.connection = await connect(
            host=self.host,
            port=self.port,
            login=self.user,
            password=self.passwd,
            virtualhost="/",
        )
        self.channel = await self.connection.channel()
        self.callback_queue = await self.channel.declare_queue(exclusive=True)
        await self.callback_queue.consume(self.on_response, no_ack=False)

        return self

    async def close(self):
        if self.connection:
            await self.connection.close()

    async def __aenter__(self):
        await self.connect()
        return self

    async def __aexit__(self, exc_type, exc_value, traceback):
        await self.close()

    async def on_response(self, message):
        if message.correlation_id is None:
            return

        callid, seqid = message.correlation_id.split(":", 1)
        seqid = int(seqid)

        if callid not in self.futures:
            return

        result = self.parse_response(message.body)
        self.results[callid][1][seqid] = result
        self.results[callid][0] -= 1

        if callid not in self.progress_bars:
            total_tasks = len(self.results[callid][1])
            self.progress_bars[callid] = tqdm(
                total=total_tasks, desc="Folding seqs", unit="requests"
            )

        self.progress_bars[callid].update(1)

        if self.results[callid][0] <= 0:
            future = self.futures.pop(callid)
            call_results = self.results.pop(callid)[1]
            future.set_result(call_results)
            self.progress_bars[callid].close()
            self.progress_bars.pop(callid, None)

    @staticmethod
    def parse_response(response):
        bp_dtype = [("i", "i4"), ("j", "i4"), ("prob", "f8")]
        free_energy, seqlen = unpack("di", response[:12])
        has_bpp = len(response) > 12 + seqlen

        if not has_bpp:
            results = {
                "structure": response[12 : 12 + seqlen].decode(),
                "free_energy": free_energy,
            }
            return results

        results = {
            "structure": response[12 : 12 + seqlen].decode(),
            "free_energy": free_energy,
            "bpp": np.frombuffer(response[12 + seqlen :], dtype=bp_dtype),
        }

        base_pair_i = results["bpp"]["i"]
        base_pair_j = results["bpp"]["j"]
        base_pair_prob = results["bpp"]["prob"]
        pair_prob_one_sided = sparse.coo_matrix(
            (base_pair_prob, (base_pair_i, base_pair_j)), shape=(seqlen, seqlen)
        )
        pair_prob = pair_prob_one_sided + pair_prob_one_sided.T
        results["pi_array"] = pair_prob

        return results

    async def call_rpc(self, payloads, loop):
        loop = loop or asyncio.get_running_loop()

        future = loop.create_future()
        callid = uuid.uuid4().hex

        self.futures[callid] = future
        self.results[callid] = [len(payloads), [None] * len(payloads)]

        # Ensure the channel is open
        if self.channel.is_closed:
            await self.connect()

        for seqid, payload in enumerate(payloads):
            message = Message(
                payload,
                reply_to=self.callback_queue.name,
                correlation_id=f"{callid}:{seqid}",
            )

            await self.channel.default_exchange.publish(message, routing_key=self.queue)

        return await future

    # NOTE: EXEUCTOR ARGUMENT EXISTS FOR COMPATIBILITY WITH LOCAL VAXIFOLD

    def call_viennarna_fold(self, seqs, executor=None):
        payloads = [
            json.dumps({"method": "viennarna_fold", "args": {"seq": seq}}).encode()
            for seq in seqs
        ]
        return self.call_rpc(payloads, None)

    def call_linearfold(self, seqs, executor=None):
        payloads = [
            json.dumps({"method": "linearfold", "args": {"seq": seq}}).encode()
            for seq in seqs
        ]
        return self.call_rpc(payloads, None)

    def call_viennarna_partition(self, seqs, executor=None):
        payloads = [
            json.dumps({"method": "viennarna_partition", "args": {"seq": seq}}).encode()
            for seq in seqs
        ]
        return self.call_rpc(payloads, None)

    def call_linearpartition(self, seqs, executor=None):
        payloads = [
            json.dumps({"method": "linearpartition", "args": {"seq": seq}}).encode()
            for seq in seqs
        ]
        return self.call_rpc(payloads, None)

    def fold(self, seqs, executor=None):
        if self.folding_engine == "viennarna":
            return self.call_viennarna_fold(seqs)
        elif self.folding_engine == "linearfold":
            return self.call_linearfold(seqs)
        else:
            raise ValueError(f"Unknown engine: {self.folding_engine}")

    def partition(self, seqs, executor=None):
        if self.partition_engine == "viennarna":
            return self.call_viennarna_partition(seqs)
        elif self.partition_engine == "linearpartition":
            return self.call_linearpartition(seqs)
        else:
            raise ValueError(f"Unknown engine: {self.partition_engine}")
