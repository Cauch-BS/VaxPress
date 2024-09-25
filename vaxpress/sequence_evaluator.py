#
# VaxPress
#
# Copyright 2023 Seoul National University
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# “Software”), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
# NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

import sys
import re
import os
import json
import uuid
import pylru
<<<<<<< HEAD
import numpy as np
from tqdm import tqdm
from concurrent import futures
from collections import Counter
from scipy import sparse
=======
import asyncio
import numpy as np
from struct import unpack
from tqdm import tqdm
from concurrent import futures
from collections import Counter
from aio_pika import Message, connect
>>>>>>> upstream/vaxifold-client
from .log import hbar_stars, log


class Evaluator():
    def __init__(self):
        self.pat_find_loops = re.compile(r'\.{2,}')

    @staticmethod
    def find_stems(structure):
        stack = []
        stemgroups = []

        for i, s in enumerate(structure):
            if s == '(':
                stack.append(i)
            elif s == ')':
                assert len(stack) >= 1
                peer = stack.pop()
                if (stemgroups and peer + 1 == stemgroups[-1][0][-1] and
                        i - 1 == stemgroups[-1][1][-1]):
                    stemgroups[-1][0].append(peer)
                    stemgroups[-1][1].append(i)
                else:
                    stemgroups.append(([peer], [i]))

        return stemgroups
    
    @staticmethod
    def unfold_unstable_structure(folding, stems):
        # TODO: This needs to be revised based on the thermodynamic model of RNA
        # folding later.
        lonepairs = [p for p in stems if len(p[0]) == 1]
        if not lonepairs:
            return folding, stems

        folding = list(folding)
        for p5, p3 in lonepairs:
            folding[p5[0]] = '.'
            folding[p3[0]] = '.'
        newstems = [p for p in stems if len(p[0]) > 1]

        return ''.join(folding), newstems

class FoldEvaluator(Evaluator):

<<<<<<< HEAD
    def __init__(self, engine: str):
        super().__init__()
        self.engine = engine
        self.initialize()
=======
    def __init__(self):
        self.pat_find_loops = re.compile(r'\.{2,}')
>>>>>>> upstream/vaxifold-client

    def __call__(self, seqs, foldings):
        for seq, folding in zip(seqs, foldings):
            # XXX: free_energy is actually the MEA, but we use it as MFE for now.
            folding['mfe'] = folding['free_energy']
            mfe_str = folding['structure']
            stems = self.find_stems(mfe_str)
            mfe_str, stems = self.unfold_unstable_structure(mfe_str, stems)
            loops = dict(Counter(map(len, self.pat_find_loops.findall(mfe_str))))

<<<<<<< HEAD
    def __call__(self, seq):
        folding, mfe = self._fold(seq)
        stems = Evaluator.find_stems(folding)
        folding, stems = Evaluator.unfold_unstable_structure(folding, stems)
        loops = dict(Counter(map(len, self.pat_find_loops.findall(folding))))
=======
            folding['folding'] = mfe_str
            folding['stems'] = stems
            folding['loops'] = loops
>>>>>>> upstream/vaxifold-client

            folding.pop('structure')
            folding.pop('free_energy')

        return foldings

class PairingProbEvaluator(Evaluator):
    def __init__(self):
        super().__init__()
        self.initialize()
    
    def initialize(self):
        try:
            import linearpartition
        except ImportError:
            raise ImportError('LinearPartition module is not available. Try "'
                              'pip install linearpartition-unofficial" to install.')
        self._pairing_prob = linearpartition.partition
    
    def __call__(self, seq):
        pred = self._pairing_prob(seq)
        bpmtx = pred['bpp']
        fe = pred['free_energy']
        mea = pred['structure']
        pi_coarray = self.get_pairingprobs(bpmtx, seq)
        stems = Evaluator.find_stems(mea)
        mea, stems = Evaluator.unfold_unstable_structure(mea, stems)
        loops = dict(Counter(map(len, self.pat_find_loops.findall(mea))))

        return {
            'pi_array': pi_coarray,
            'efe': fe, #efe = ensemble free energy
            'folding': mea, #mea = maximum expected accuracy
            'stems': stems,
            'loops': loops,
        }
    
    @staticmethod
    def get_pairingprobs(bpmtx, seq):
        # returns scipy sparse matrix (coo_array)
        I = bpmtx['i']
        J = bpmtx['j']
        prob = bpmtx['prob']
        A = sparse.coo_array((prob, (I, J)), shape=(len(seq), len(seq)))
        return (A+A.T)

class SequenceEvaluator:

    folding_cache_size = 8192
    bpp_cache_size = 8192

    def __init__(self, scoring_funcs, scoreopts, execopts, mutantgen, species,
                 length_cds, quiet):
        self.scoring_funcs = scoring_funcs
        self.scoreopts = scoreopts
        self.execopts = execopts

        self.length_cds = length_cds
        self.mutantgen = mutantgen
        self.species = species

        self.quiet = quiet

        self.initialize()

    def initialize_vaxifold_client(self):
        vaxifold_url = os.environ['VAXIFOLD_URL']
        vaxifold_queue = os.environ['VAXIFOLD_QUEUE']
        engine = self.execopts.folding_engine

        self.vaxifold = AsyncVaxiFoldClient(vaxifold_url, vaxifold_queue)
        if engine == 'viennarna':
            self._partition = self.vaxifold.call_viennarna_partition
        elif engine == 'linearpartition':
            self._partition = self.vaxifold.call_linearpartition
        else:
            raise ValueError(f'Unsupported RNA folding engine: {engine}')

        asyncio.run(self.vaxifold.connect())

    def initialize(self):
<<<<<<< HEAD
        self.foldeval = FoldEvaluator(self.execopts.folding_engine)
        self.bpp_eval = PairingProbEvaluator()
=======
        self.initialize_vaxifold_client()

        self.foldeval = FoldEvaluator()
>>>>>>> upstream/vaxifold-client
        self.folding_cache = pylru.lrucache(self.folding_cache_size)
        self.bpp_cache = pylru.lrucache(self.bpp_cache_size)

        self.scorefuncs_nofolding = []
        self.scorefuncs_folding = []
        self.scorefuncs_bpp =[]
        self.annotationfuncs = []
        self.penalty_metric_flags = {}

        additional_opts = {
            '_length_cds': self.length_cds,
        }

        use_fold, use_mea = True, True
        for funcname, cls in self.scoring_funcs.items():
            opts = self.scoreopts[funcname]
            if 'mfe' in funcname:
                if ('weight' in opts and opts['weight'] == 0):
                    use_fold = False
            if 'aup' in funcname:
                if ('weight' in opts and opts['weight'] == 0):
                    use_mea = False

        for funcname, cls in self.scoring_funcs.items():
            funcoff = False
            opts = self.scoreopts[funcname]
            if (('weight' in opts and opts['weight'] == 0) or
                    ('off' in opts and opts['off'])):
                if not cls.use_annotation_on_zero_weight:
                    continue
                funcoff = True

            opts.update(additional_opts)
            for reqattr in cls.requires:
                opts['_' + reqattr] = getattr(self, reqattr)

            try:
                scorefunc_inst = cls(**opts)
            except EOFError:
                continue

            self.annotationfuncs.append(scorefunc_inst)
            if funcoff:
                continue

            if cls.uses_folding and use_fold:
                self.scorefuncs_folding.append(scorefunc_inst)
            elif cls.uses_basepairing_prob and use_mea and (
                not cls.uses_folding or not use_fold):
                self.scorefuncs_bpp.append(scorefunc_inst)       
            else:
                self.scorefuncs_nofolding.append(scorefunc_inst)

            self.penalty_metric_flags.update(cls.penalty_metric_flags)
        # log.info('>>>>Scoring functions used:')
        # log.info('\t Uses Folding: {}'.format(self.scorefuncs_folding))
        # log.info('\t Uses Base Pairing Probability: {}'.format(self.scorefuncs_bpp))
        # log.info('\t Does not use Folding: {}'.format(self.scorefuncs_nofolding))

    def evaluate(self, seqs, executor):
        with SequenceEvaluationSession(self, seqs, executor) as sess:
            sess.evaluate()

            if not sess.errors:
                total_scores = [sum(s.values()) for s in sess.scores]
                return total_scores, sess.scores, sess.metrics, sess.foldings, sess.pairingprobs
            else:
                return None, None, None, None, None

    async def fold_individual_seq(self, seq):
        foldings = await self._partition([{seq: seq}])
        annotated = self.foldeval([seq], foldings)
        self.folding_cache[seq] = annotated[0]

    def get_folding(self, seq):
        if seq not in self.folding_cache:
            asyncio.run(self.fold_individual_seq(seq))
        return self.folding_cache[seq]
    
    def get_pairing_prob(self, seq):
        if seq not in self.bpp_cache:
            self.bpp_cache[seq] = self.bpp_eval(seq)
        return self.bpp_cache[seq]

    def prepare_evaluation_data(self, seq):
        folding = self.get_folding(seq)
        pairingprob = self.get_pairing_prob(seq)

        seqevals = {}
        seqevals['local-metrics'] = localmet = {}
        for fun in self.annotationfuncs:
            if hasattr(fun, 'evaluate_local'):
                if fun.uses_folding:
                    localmet.update(fun.evaluate_local(seq, folding))
                elif fun.uses_basepairing_prob:
                    localmet.update(fun.evaluate_local(seq, pairingprob))
                else:
                    localmet.update(fun.evaluate_local(seq))

            if hasattr(fun, 'annotate_sequence'):
                if fun.uses_folding:
                    seqevals.update(fun.annotate_sequence(seq, folding))
                elif fun.uses_basepairing_prob:
                    seqevals.update(fun.annotate_sequence(seq, pairingprob))
                else:
                    seqevals.update(fun.annotate_sequence(seq))

        return seqevals


class AsyncVaxiFoldClient:

    def __init__(self, url, queue):
        self.url = url
        self.queue = queue

        self.futures = {}
        self.results = {}

    async def connect(self):
        self.connection = await connect(self.url)
        self.channel = await self.connection.channel()
        self.callback_queue = await self.channel.declare_queue(exclusive=True)
        await self.callback_queue.consume(self.on_response, no_ack=False)

        return self

    async def on_response(self, message):
        if message.correlation_id is None:
            return

        callid, seqid = message.correlation_id.split(':', 1)
        seqid = int(seqid)

        if callid not in self.futures:
            return

        result = self.parse_partition_response(message.body)
        self.results[callid][1][seqid] = result
        self.results[callid][0] -= 1
        print('.', end='', flush=True)

        if self.results[callid][0] <= 0:
            future = self.futures.pop(callid)
            call_results = self.results.pop(callid)[1]
            future.set_result(call_results)
            print(' done')

    @staticmethod
    def parse_partition_response(response):
        bp_dtype = [('i', 'i4'), ('j', 'i4'), ('prob', 'f8')]
        free_energy, seqlen = unpack('di', response[:12])
        return {
            'structure': response[12:12+seqlen].decode(),
            'free_energy': free_energy,
            'bpp': np.frombuffer(response[12+seqlen:], dtype=bp_dtype),
        }

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
                payload, reply_to=self.callback_queue.name,
                correlation_id=f'{callid}:{seqid}'
            )

            await self.channel.default_exchange.publish(
                message, routing_key=self.queue)

        return await future

    def call_viennarna_partition(self, seqs, loop=None):
        payloads = [
            json.dumps({'method': 'viennarna_partition', 'args': {'seq': seq}}).encode()
            for seq in seqs]
        return self.call_rpc(payloads, loop)

    def call_linearpartition(self, seqs, loop=None):
        payloads = [
            json.dumps({'method': 'linearpartition', 'args': {'seq': seq}}).encode()
            for seq in seqs]
        return self.call_rpc(payloads, loop)


class SequenceEvaluationSession:

    def __init__(self, evaluator: SequenceEvaluator, seqs: list[str],
                 executor: futures.Executor):
        self.seqs = seqs
        self.executor = executor

        self.scores = [{} for i in range(len(seqs))]
        self.metrics = [{} for i in range(len(seqs))]
        self.foldings = [None] * len(seqs)
        self.pairingprobs = [None] * len(seqs)
        self.errors = []

        self.foldeval = evaluator.foldeval
        self.vaxifold_partition = evaluator._partition

        self.bpp_cache = evaluator.bpp_cache
        self.bpps_remaining = len(seqs)
        self.bpp_eval = evaluator.bpp_eval

        self.num_tasks = (
            len(evaluator.scorefuncs_folding) +
            len(evaluator.scorefuncs_nofolding) +
            len(evaluator.scorefuncs_bpp) + 
            len(seqs))

        self.quiet = evaluator.quiet

        self.scorefuncs_folding = evaluator.scorefuncs_folding
        self.scorefuncs_nofolding = evaluator.scorefuncs_nofolding
        self.scorefuncs_bpp = evaluator.scorefuncs_bpp
        self.annotationfuncs = evaluator.annotationfuncs

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def evaluate(self):
        asyncio.run(self.run_evaluations())

<<<<<<< HEAD
        # First, base pairing probability is calculated.
        for i, seq in enumerate(self.seqs):
            if self.errors:
                continue

            if seq in self.bpp_cache:
                self.pairingprobs[i] = self.bpp_cache[seq]
                self.bpps_remaining -= 1
                if self.pbar is not None:
                    self.pbar.update()
                continue

            future = self.executor.submit(self.bpp_eval, seq)
            future._seqidx = i
            future._type = 'partition'
            jobs.add(future)

        # Secondary structure prediction is the second set of tasks.
        for i, seq in enumerate(self.seqs):
            if self.errors: # skip remaining tasks on error
                continue
=======
    async def run_evaluations(self):
        tasks_with_folding = self.run_scoring_with_folding()
        tasks_without_folding = self.run_scoring_without_folding()
>>>>>>> upstream/vaxifold-client

        await asyncio.gather(tasks_without_folding, tasks_with_folding)

<<<<<<< HEAD
            future = self.executor.submit(self.foldeval, seq)
            future._seqidx = i
            future._type = 'folding'
            jobs.add(future)

        # Then, scoring functions that does not require folding are executed.
        for scorefunc in self.scorefuncs_nofolding: # no folding, no bpp
            if self.errors:
                continue
=======
    async def run_scoring_without_folding(self):
        tasks = []
        for scorefunc in self.scorefuncs_nofolding:
            task = self.call_scorefunc(scorefunc, self.seqs)
            tasks.append(task)
>>>>>>> upstream/vaxifold-client

        await asyncio.gather(*tasks)

<<<<<<< HEAD
        # Wait until all folding tasks are finished.
        while jobs and not self.errors and (self.foldings_remaining > 0
                                            or self.bpps_remaining > 0):
            done, jobs = futures.wait(jobs, timeout=0.1,
                                      return_when=futures.FIRST_COMPLETED)
            for future in done:
                if future._type == 'folding':
                    self.collect_folding(future)
                elif future._type == 'partition':
                    self.collect_pairingprob(future)
                elif future._type == 'scoring':
                    self.collect_scores(future)
        
        # Scoring functions requiring base pairing probability are executed.
        for scorefunc in self.scorefuncs_bpp:
            if self.errors:
                continue

            future = self.executor.submit(scorefunc, self.seqs,
                                          pairingprobs = self.pairingprobs)
            future._type = 'scoring'
            jobs.add(future)

        # Scoring functions requiring mfe are executed.
        for scorefunc in self.scorefuncs_folding:
            if self.errors:
                continue

            future = self.executor.submit(scorefunc, self.seqs,
                                          foldings = self.foldings)
            future._type = 'scoring'
            jobs.add(future)
        
        while jobs and not self.errors:
            done, jobs = futures.wait(jobs, timeout=0.1)
            for future in done:
                if future._type == 'folding':
                    self.collect_folding(future)
                elif future._type == 'partition':
                    self.collect_pairingprob(future)
                elif future._type == 'scoring':
                    self.collect_scores(future)

    def collect_scores(self, future):
        try:
            ret = future.result()
            #print(ret)
=======
    async def run_scoring_with_folding(self):
        loop = asyncio.get_running_loop()

        foldings = await self.vaxifold_partition(self.seqs)
        future = loop.run_in_executor(None, self.foldeval, self.seqs, foldings)
        try:
            self.foldings = await future
        except Exception as exc:
            return self.handle_exception(exc)

        tasks = []
        for scorefunc in self.scorefuncs_folding:
            task = self.call_scorefunc(scorefunc, self.seqs, self.foldings)
            tasks.append(task)

        await asyncio.gather(*tasks)

    async def call_scorefunc(self, scorefunc, seqs, foldings=None):
        loop = asyncio.get_running_loop()

        if foldings is None:
            future = loop.run_in_executor(None, scorefunc, seqs)
        else:
            future = loop.run_in_executor(None, scorefunc, seqs, foldings)

        try:
            ret = await future
>>>>>>> upstream/vaxifold-client
            if ret is None:
                self.errors.append('KeyboardInterrupt')
                return
            scoreupdates, metricupdates = ret
        except Exception as exc:
            return self.handle_exception(exc)

        # Update scores
        for k, updates in scoreupdates.items():
            assert len(updates) == len(self.scores)
            for s, u in zip(self.scores, updates):
                s[k] = u

        # Update metrics
        for k, updates in metricupdates.items():
            assert len(updates) == len(self.metrics)
            for s, u in zip(self.metrics, updates):
                s[k] = u
    
    def collect_pairingprob(self, future):
        try:
            pairingprob = future.result()
            if pairingprob is None:
                self.errors.append('KeyboardInterrupt')
                if self.pbar is not None:
                    self.pbar.close()
                self.pbar = None
                return
        except Exception as exc:
            return self.handle_exception(exc)
        i = future._seqidx
        self.pairingprobs[i] = pairingprob
        self.bpp_cache[self.seqs[i]] = pairingprob
        self.bpps_remaining -= 1

        if self.pbar is not None:
            self.pbar.update()

    def handle_exception(self, exc):
        import traceback
        import io

        errormsg = io.StringIO()
        traceback.print_exc(file=errormsg)

        msg = [
            hbar_stars,
            'Error occurred in a scoring function:',
            errormsg.getvalue(),
            hbar_stars,
            '',
            'Termination in progress. Waiting for running tasks '
            'to finish before closing the program.']

        log.error('\n'.join(msg))
        self.errors.append(exc.args)
