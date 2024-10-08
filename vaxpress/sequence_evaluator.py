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
import os
import re
import pylru  # type: ignore
import asyncio
from concurrent import futures
from collections import Counter
from .vaxifold import AsyncVaxiFoldClient
from .log import hbar_stars, log
from tqdm import tqdm


class ParseStructure:

    def __init__(self):
        self.pat_find_loops = re.compile(r"\.{2,}")

    @staticmethod
    def find_stems(structure):
        stack = []
        stemgroups = []

        for i, s in enumerate(structure):
            if s == "(":
                stack.append(i)
            elif s == ")":
                assert len(stack) >= 1
                peer = stack.pop()
                if (
                    stemgroups
                    and peer + 1 == stemgroups[-1][0][-1]
                    and i - 1 == stemgroups[-1][1][-1]
                ):
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
            folding[p5[0]] = "."
            folding[p3[0]] = "."
        newstems = [p for p in stems if len(p[0]) > 1]

        return "".join(folding), newstems

    def __call__(self, seqs, foldings):
        for _, folding in zip(seqs, foldings):
            # XXX: free_energy is actually the MEA, but we use it as MFE for now.
            folding["mfe"] = folding["free_energy"]
            mfe_str = folding["structure"]
            stems = self.find_stems(mfe_str)
            mfe_str, stems = self.unfold_unstable_structure(mfe_str, stems)
            loops = dict(Counter(map(len, self.pat_find_loops.findall(mfe_str))))

            folding["folding"] = mfe_str
            folding["stems"] = stems
            folding["loops"] = loops

            folding.pop("free_energy")
            folding.pop("structure")

        return foldings


class SequenceEvaluator:

    folding_cache_size = 16384

    def __init__(
        self,
        scoring_funcs,
        scoreopts,
        execopts,
        mutantgen,
        species,
        length_cds,
        quiet,
        host=None,
        port=None,
        user=None,
        passwd=None,
        queue=None,
    ):
        self.scoring_funcs = scoring_funcs
        self.scoreopts = scoreopts
        self.execopts = execopts

        self.length_cds = length_cds
        self.mutantgen = mutantgen
        self.species = species

        self.quiet = quiet
        self.use_fold = True
        self.foldeval = ParseStructure()
        self.folding_cache = pylru.lrucache(self.folding_cache_size)

        self.scorefuncs_nofolding = []
        self.scorefuncs_folding = []
        self.scorefuncs_bpp = []
        self.annotationfuncs = []
        self.penalty_metric_flags = {}
        self.host = host
        self.port = port
        self.user = user
        self.passwd = passwd
        self.queue = queue

        self.initialize()

    def initialize(self):

        host = self.host
        port = self.port
        user = self.user
        passwd = self.passwd
        queue = self.queue

        if not host:
            if "RABBITMQ_HOST" in os.environ:
                host = os.environ["RABBITMQ_HOST"]
            else:
                raise ConnectionError("RabbitMQ host not specified.")
        if not port:
            if "RABBITMQ_PORT" in os.environ:
                port = os.environ["RABBITMQ_PORT"]
            else:
                raise ConnectionError("RabbitMQ port not specified.")
        if not user:
            if "RABBITMQ_USER" in os.environ:
                user = os.environ["RABBITMQ_USER"]
            else:
                raise ConnectionError("RabbitMQ user not specified.")
        if not passwd:
            if "RABBITMQ_PASS" in os.environ:
                passwd = os.environ["RABBITMQ_PASS"]
            else:
                raise ConnectionError("RabbitMQ password not specified.")
        if queue is None:
            if "VAXIFOLD_QUEUE" in os.environ:
                queue = os.environ["VAXIFOLD_QUEUE"]
            else:
                raise ConnectionError("VaxiFold queue not specified.")

        self.vaxifold = AsyncVaxiFoldClient(
            host=host,
            port=port,
            user=user,
            passwd=passwd,
            queue=queue,
            folding_engine=self.execopts.folding_engine,
            partition_engine=self.execopts.partition_engine,
        )

        additional_opts = {
            "_length_cds": self.length_cds,
        }

        for funcname, cls in self.scoring_funcs.items():
            opts = self.scoreopts[funcname]
            if "mfe" in funcname:
                if "weight" in opts and opts["weight"] == 0:
                    self.use_fold = False

        for funcname, cls in self.scoring_funcs.items():
            funcoff = False
            opts = self.scoreopts[funcname]
            if ("weight" in opts and opts["weight"] == 0) or (
                "off" in opts and opts["off"]
            ):
                if not cls.use_annotation_on_zero_weight:
                    continue
                funcoff = True

            opts.update(additional_opts)
            for reqattr in cls.requires:
                opts["_" + reqattr] = getattr(self, reqattr)

            try:
                scorefunc_inst = cls(**opts)
            except EOFError:
                continue

            self.annotationfuncs.append(scorefunc_inst)
            if funcoff:
                continue

            if cls.uses_folding and self.use_fold:
                self.scorefuncs_folding.append(scorefunc_inst)
            elif cls.uses_basepairing_prob or (cls.uses_folding and not self.use_fold):
                self.scorefuncs_bpp.append(scorefunc_inst)
            else:
                self.scorefuncs_nofolding.append(scorefunc_inst)

            self.penalty_metric_flags.update(cls.penalty_metric_flags)

    def find_fitness(self, seqs, executor):
        with FitnessEvaluationSession(self, seqs, executor) as f_sess:
            f_sess.evaluate()

            if not f_sess.errors:
                total_scores = [sum(s.values()) for s in f_sess.scores]
                return (
                    total_scores,
                    f_sess.scores,
                    f_sess.metrics,
                    f_sess.foldings,
                )
            else:
                return None, None, None, None

    async def parse_individual_seq(self, seq):
        try:
            async with self.vaxifold:
                if self.use_fold:
                    mfe_foldings = await self.vaxifold.fold([seq])
                pf_foldings = await self.vaxifold.partition([seq])
        except Exception as exc:
            log.error(f"Error occurred in folding: {exc}")

        pf_annotated = self.foldeval([seq], pf_foldings)
        if self.use_fold:
            mfe_annotated = self.foldeval([seq], mfe_foldings)
        else:
            mfe_annotated = [dict()]

        self.folding_cache[seq] = mfe_annotated[0], pf_annotated[0]

    def parse_seq(self, seq):
        if seq not in self.folding_cache:
            asyncio.run(self.parse_individual_seq(seq))
        return self.folding_cache[seq]

    def prepare_evaluation_data(self, seq):
        mfe_parsed, pf_parsed = self.parse_seq(seq)
        seqevals = {}
        seqevals["local-metrics"] = localmet = {}
        for fun in self.annotationfuncs:
            if hasattr(fun, "evaluate_local"):
                if fun.uses_folding and self.use_fold:
                    localmet.update(fun.evaluate_local(seq, mfe_parsed))
                elif fun.uses_basepairing_prob or (
                    fun.uses_folding and not self.use_fold
                ):
                    localmet.update(fun.evaluate_local(seq, pf_parsed))
                else:
                    localmet.update(fun.evaluate_local(seq))

            if hasattr(fun, "annotate_sequence"):
                if fun.uses_folding and self.use_fold:
                    seqevals.update(fun.annotate_sequence(seq, mfe_parsed))
                elif fun.uses_basepairing_prob or (
                    fun.uses_folding and not self.use_fold
                ):
                    seqevals.update(fun.annotate_sequence(seq, pf_parsed))
                else:
                    seqevals.update(fun.annotate_sequence(seq))

        return seqevals


class FitnessEvaluationSession:

    def __init__(
        self, evaluator: SequenceEvaluator, seqs: list[str], executor: futures.Executor
    ):
        self.seqs = seqs
        self.executor = executor

        self.scores: list[dict] = [{} for _ in range(len(seqs))]
        self.metrics: list[dict] = [{} for _ in range(len(seqs))]
        self.foldings = [None] * len(seqs)
        self.errors: list = []

        self.foldeval = evaluator.foldeval
        self.folding_cache = evaluator.folding_cache

        self.vaxifold = evaluator.vaxifold

        self.pbar_nofold = None
        self.pbar_fold = None
        self.quiet = evaluator.quiet

        self.scorefuncs_folding = evaluator.scorefuncs_folding
        self.scorefuncs_nofolding = evaluator.scorefuncs_nofolding
        self.scorefuncs_bpp = evaluator.scorefuncs_bpp
        self.annotationfuncs = evaluator.annotationfuncs

        self.scorefuncs_vaxifold = self.scorefuncs_folding[:]

        for item in self.scorefuncs_bpp:
            if item not in self.scorefuncs_vaxifold:
                self.scorefuncs_vaxifold.append(item)

        self.num_tasks_fold = len(self.scorefuncs_vaxifold)
        self.num_tasks_nofold = len(self.scorefuncs_nofolding)

    def __enter__(self):
        log.info("")
        self.pbar_nofold = tqdm(
            total=self.num_tasks_nofold,
            desc="Fitness Eval (No Fold)",
            file=sys.stderr,
            unit="task",
            disable=self.quiet,
        )
        self.pbar_fold = tqdm(
            total=self.num_tasks_fold,
            desc="Fitness Eval (Fold)",
            file=sys.stderr,
            unit="task",
            disable=self.quiet,
        )
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.pbar_nofold:
            self.pbar_nofold.close()
        if self.pbar_fold:
            self.pbar_fold.close()
        pass

    async def call_scorefunc(self, scorefunc, seqs, foldings=None):
        loop = asyncio.get_running_loop()

        if foldings is None:
            future = loop.run_in_executor(None, scorefunc, seqs)
        else:
            future = loop.run_in_executor(None, scorefunc, seqs, foldings)

        try:
            ret = await future
            if ret is None:
                self.errors.append("KeyboardInterrupt")
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

    async def run_scoring_without_folding(self):
        tasks = []
        for scorefunc in self.scorefuncs_nofolding:
            task = self.call_scorefunc(scorefunc, self.seqs)
            tasks.append(task)
            self.pbar_nofold.update(1)
        await asyncio.gather(*tasks)

    async def run_scoring_with_folding(self):
        loop = asyncio.get_running_loop()

        mfe_results = [None for _ in range(len(self.seqs))]
        pf_results = [None for _ in range(len(self.seqs))]

        new_seqs = []

        for i, seq in enumerate(self.seqs):
            if seq in self.folding_cache:
                mfe_results[i], pf_results[i] = self.folding_cache[seq]
                continue
            new_seqs.append(seq)

        try:
            async with self.vaxifold:
                if self.scorefuncs_folding:
                    mfe_foldings = await self.vaxifold.fold(new_seqs)
                pf_foldings = await self.vaxifold.partition(new_seqs)

        except Exception as exc:
            return self.handle_exception(exc)

        pf_future = loop.run_in_executor(None, self.foldeval, self.seqs, pf_foldings)

        if self.scorefuncs_folding:
            mfe_future = loop.run_in_executor(
                None, self.foldeval, self.seqs, mfe_foldings
            )
            try:
                new_mfe_results = await mfe_future
                new_pf_results = await pf_future
                for i, entry in enumerate(self, zip(mfe_results, pf_results)):
                    mfe_entry, pf_entry = entry
                    if (mfe_entry) and (pf_entry):
                        continue
                    mfe_results[i] = new_mfe_results.pop(0)
                    pf_results[i] = new_pf_results.pop(0)
                    if not isinstance(mfe_results[i], dict):
                        mfe_results[i] = mfe_results[i][0]
                    if not isinstance(pf_results[i], dict):
                        pf_results[i] = pf_results[i][0]
                    mfe_results[i].update(
                        {
                            key: pf_results[i][key]
                            for key in pf_results[i]
                            if key not in mfe_results[i]
                        }
                    )
                    self.folding_cache[self.seqs[i]] = mfe_results[i], pf_results[i]
                self.foldings = mfe_results

            except Exception as exc:
                return self.handle_exception(exc)

        elif not self.scorefuncs_folding:
            try:
                new_pf_results = await pf_future
                for i, pf_entry in enumerate(pf_results):
                    if pf_entry:
                        continue
                    pf_results[i] = new_pf_results.pop(0)
                    if not isinstance(pf_results[i], dict):
                        pf_results[i] = pf_results[i][0]
                    self.folding_cache[self.seqs[i]] = dict(), pf_results[i]
                self.foldings = pf_results
            except Exception as exc:
                return self.handle_exception(exc)

        tasks = []

        for scorefunc in self.scorefuncs_vaxifold:
            task = self.call_scorefunc(scorefunc, self.seqs, self.foldings)
            tasks.append(task)
            self.pbar_fold.update(1)

        await asyncio.gather(*tasks)

    async def run_evaluations(self):
        tasks_without_folding = self.run_scoring_without_folding()
        tasks_with_folding = self.run_scoring_with_folding()

        await asyncio.gather(tasks_without_folding, tasks_with_folding)

    def evaluate(self):
        asyncio.run(self.run_evaluations())

    def handle_exception(self, exc):
        import traceback
        import io

        errormsg = io.StringIO()
        traceback.print_exc(file=errormsg)

        msg = [
            hbar_stars,
            "Error occurred in a scoring function:",
            errormsg.getvalue(),
            hbar_stars,
            "",
            "Termination in progress. Waiting for running tasks "
            "to finish before closing the program.",
        ]

        log.error("\n".join(msg))
        self.errors.append(exc.args)
