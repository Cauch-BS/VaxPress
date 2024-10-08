from aio_pika import connect, Message
import json
import uuid
import asyncio
from struct import unpack
import numpy as np
from scipy import sparse  # type: ignore
from tqdm import tqdm


class AsyncVaxiFoldClient:

    def __init__(self, url, queue, folding_engine, partition_engine):
        self.url = url
        self.queue = queue
        self.folding_engine = folding_engine
        self.partition_engine = partition_engine
        self.progress_bars = {}
        self.futures = {}
        self.results = {}

    async def connect(self):
        self.connection = await connect(self.url)
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

        result = self.parse_partition_response(message.body)
        self.results[callid][1][seqid] = result
        self.results[callid][0] -= 1

        if callid not in self.progress_bars:
            total_tasks = len(self.results[callid][1])
            self.progress_bars[callid] = tqdm(
                total=total_tasks, desc=f"Folding seqs", unit="requests"
            )

        self.progress_bars[callid].update(1)

        if self.results[callid][0] <= 0:
            future = self.futures.pop(callid)
            call_results = self.results.pop(callid)[1]
            future.set_result(call_results)
            self.progress_bars[callid].close()
            self.progress_bars.pop(callid, None)

    @staticmethod
    def parse_partition_response(response):
        bp_dtype = [("i", "i4"), ("j", "i4"), ("prob", "f8")]
        free_energy, seqlen = unpack("di", response[:12])

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
        loop = loop if loop is not None else asyncio.get_running_loop()

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

    def call_viennarna_fold(self, seqs, loop=None):
        payloads = [
            json.dumps({"method": "viennarna_fold", "args": {"seq": seq}}).encode()
            for seq in seqs
        ]
        return self.call_rpc(payloads, loop)

    def call_linearfold(self, seqs, loop=None):
        payloads = [
            json.dumps({"method": "linearfold", "args": {"seq": seq}}).encode()
            for seq in seqs
        ]
        return self.call_rpc(payloads, loop)

    def call_viennarna_partition(self, seqs, loop=None):
        payloads = [
            json.dumps({"method": "viennarna_partition", "args": {"seq": seq}}).encode()
            for seq in seqs
        ]
        return self.call_rpc(payloads, loop)

    def call_linearpartition(self, seqs, loop=None):
        payloads = [
            json.dumps({"method": "linearpartition", "args": {"seq": seq}}).encode()
            for seq in seqs
        ]
        return self.call_rpc(payloads, loop)

    def fold(self, seqs, loop=None):
        if self.folding_engine == "viennarna":
            return self.call_viennarna_fold(seqs, loop)
        elif self.folding_engine == "linearfold":
            return self.call_linearfold(seqs, loop)
        else:
            raise ValueError(f"Unknown engine: {self.folding_engine}")

    def partition(self, seqs, loop=None):
        if self.partition_engine == "viennarna":
            return self.call_viennarna_partition(seqs, loop)
        elif self.partition_engine == "linearpartition":
            return self.call_linearpartition(seqs, loop)
        else:
            raise ValueError(f"Unknown engine: {self.partition_engine}")
