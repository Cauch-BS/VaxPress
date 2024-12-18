from aio_pika import connect, Message
import json
import uuid
import asyncio
from struct import unpack
import numpy as np
from tqdm.asyncio import tqdm_asyncio as atqdm  # type: ignore

from concurrent.futures import Executor, ProcessPoolExecutor
from importlib.metadata import version
from .data.thermodynamics import M1PSI_PARAMS


class LocalVaxiFold:

    def __init__(self, folding_engine: str, partition_engine: str):
        self.folding_engine = folding_engine
        self.partition_engine = partition_engine

    async def __aenter__(self):
        return self

    async def __aexit__(self, exc_type, exc_value, traceback):
        pass

    async def fold(
        self,
        seqs,
        mod: str = None,  # type: ignore[assignment]
        executor: Executor = None,  # type: ignore[assignment]
    ):
        executor = executor or ProcessPoolExecutor()
        if self.folding_engine == "viennarna":
            func = self.viennarna_fold
        elif self.folding_engine == "linearfold":
            func = self.linearfold
        else:
            raise ValueError(f"Unknown engine: {self.folding_engine}")

        total_tasks = len(seqs)
        progress_bar = atqdm(total=total_tasks, desc="Folding seqs", unit="requests")

        async def fold_with_progress(seq, mod):
            result = await asyncio.wrap_future(executor.submit(func, seq, mod))
            progress_bar.update(1)
            return result

        tasks = [fold_with_progress(seq, mod) for seq in seqs]
        results = await asyncio.gather(*tasks)

        progress_bar.close()
        return results

    async def partition(
        self,
        seqs,
        mod: str = None,  # type: ignore[assignment]
        executor: Executor = None,  # type: ignore[assignment]
    ):
        executor = executor or ProcessPoolExecutor()
        if self.partition_engine == "viennarna":
            func = self.viennarna_partition
        elif self.partition_engine == "linearpartition":
            func = self.linearpartition
        else:
            raise ValueError(f"Unknown engine: {self.partition_engine}")

        total_tasks = len(seqs)
        progress_bar = atqdm(total=total_tasks, desc="Partition seqs", unit="requests")

        async def partition_with_progress(seq, mod):
            result = await asyncio.wrap_future(executor.submit(func, seq, mod))
            progress_bar.update(1)
            return result

        tasks = [partition_with_progress(seq, mod) for seq in seqs]
        results = await asyncio.gather(*tasks)

        progress_bar.close()
        return results

    @staticmethod
    def viennarna_fold(seq, mod=None):
        try:
            import ViennaRNA as RNA  # type: ignore
        except ImportError:
            raise ImportError(
                "ViennaRNA is not installed. \n try: `pip install ViennaRNA`"
            )
        assert version("ViennaRNA") >= "2.4.0", "ViennaRNA version >= 2.4.0 is required"
        if not mod:
            struct, mfe = RNA.fold(seq)
            return {"structure": struct, "free_energy": mfe}

        fc = RNA.fold_compound(seq)
        fc.sc_mod(
            params=M1PSI_PARAMS,
            modification_sites=[i for i in range(1, len(seq) + 1) if seq[i - 1] == "U"],
        )
        struct, mfe = fc.mfe()
        return {"structure": struct, "free_energy": mfe}

    @staticmethod
    def linearfold(seq, mod=None):
        try:
            import linearfold  # type: ignore
        except ImportError:
            raise ImportError(
                "LinearFold is not installed. \n try: `pip install linearfold-unofficial`"
            )
        assert (
            version("linearfold-unofficial") >= "0.1"
        ), "LinearFold version >= 0.1 is required"
        if mod:
            raise NotImplementedError("LinearFold does not support modifications")
        struct, mfe = linearfold.fold(seq)
        return {"structure": struct, "free_energy": mfe}

    @staticmethod
    def viennarna_partition(seq, mod=None):
        try:
            import ViennaRNA as RNA  # type: ignore
        except ImportError:
            raise ImportError(
                "ViennaRNA is not installed. \n try: `pip install ViennaRNA`"
            )
        assert version("ViennaRNA") >= "2.4.0", "ViennaRNA version >= 2.4.0 is required"

        fc = RNA.fold_compound(seq)
        if mod:
            fc.sc_mod(
                params=M1PSI_PARAMS,
                modification_sites=[
                    i for i in range(1, len(seq) + 1) if seq[i - 1] == "U"
                ],
            )
        _, fe = fc.pf()
        fold, _ = fc.MEA()
        bpp = np.array(fc.bpp())[1:, 1:]
        pi_array = np.sum(bpp + bpp.T, axis=0)
        del bpp
        results = {
            "structure": fold,
            "free_energy": fe,
            "pi_array": pi_array,
        }
        assert len(fold) == len(seq), "Partitioning failed for sequence: %s" % seq
        return results

    @staticmethod
    def linearpartition(seq, mod=None):
        try:
            import linearpartition  # type: ignore
        except ImportError:
            raise ImportError(
                "LinearPartition is not installed. \n try: `pip install linearpartition-unofficial`"
            )
        assert (
            version("linearpartition-unofficial") >= "0.3"
        ), "LinearPartition version >= 0.3 is required"
        if mod:
            raise NotImplementedError("LinearPartition does not support modifications")
        pred = linearpartition.partition(seq)
        results = {
            "structure": pred["structure"],
            "free_energy": pred["free_energy"],
        }
        probbypos = np.zeros(len(seq), dtype=np.float64)
        np.add.at(probbypos, pred["bpp"]["i"], pred["bpp"]["prob"])
        np.add.at(probbypos, pred["bpp"]["j"], pred["bpp"]["prob"])
        pi_array = probbypos.clip(0, 1)
        results["pi_array"] = pi_array

        del pred
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
            self.progress_bars[callid] = atqdm(
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
        }
        bpp = np.frombuffer(response[12 + seqlen :], dtype=bp_dtype)
        probbypos = np.zeros(seqlen, dtype=np.float64)
        probbypos[bpp["i"]] = bpp["prob"]
        probbypos[bpp["j"]] += bpp["prob"]
        pi_array = probbypos.clip(0, 1)

        results["pi_array"] = pi_array

        del bpp

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
