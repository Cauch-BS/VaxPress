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

import os
import sys
import time
from collections import namedtuple
from concurrent import futures
from itertools import cycle

# import vaxpress.nsga2
from textwrap import wrap
from typing import Any, Dict, Iterator

import numpy as np
from tabulate import tabulate  # type: ignore[import-untyped]

from . import __version__
from .log import hbar, hbar_double, log
from .mutant_generator import MutantGenerator
from .presets import dump_to_preset
from .sequence import Sequence
from .sequence_evaluator import SequenceEvaluator

ExecutionOptions = namedtuple(
    "ExecutionOptions",
    [
        "n_iterations",
        "n_population",
        "n_survivors",
        "initial_mutation_rate",
        "winddown_trigger",
        "winddown_rate",
        "output",
        "command_line",
        "overwrite",
        "seed",
        "processes",
        "random_initialization",
        "conservative_start",
        "boost_loop_mutations",
        "full_scan_interval",
        "species",
        "codon_table",
        "protein",
        "is_modified",
        "preserve_stop",
        "quiet",
        "seq_description",
        "print_top_mutants",
        "addons",
        "lineardesign_dir",
        "lineardesign_lambda",
        "lineardesign_omit_start",
        "lineardesign_penalty",
        "lineardesign_penalty_weight",
        "folding_engine",
        "partition_engine",
        "host",
        "port",
        "user",
        "passwd",
        "queue",
    ],
)


class CDSEvolutionChamber:

    stop_threshold = 0.2
    table_header_length = 6

    finalized_seqs = None
    checkpoint_path = None
    fasta_line_width = 72

    def __init__(
        self,
        rawseqs: list[str],
        cdsseqs: list[str],
        utr5s: list[str],
        utr3s: list[str],
        scoring_funcs: dict,
        scoring_options: dict,
        exec_options: ExecutionOptions,
    ):
        self.seqall = [
            Sequence(
                seq,
                exec_options.codon_table,
                exec_options.protein,
                cdsseq,
                utr5,
                utr3,
                exec_options.preserve_stop,
            )
            for seq, cdsseq, utr5, utr3 in zip(rawseqs, cdsseqs, utr5s, utr3s)
        ]
        assert len(set((seq.translate(seq.cdsseq) for seq in self.seqall))) == 1, (
            "The sequences must have the same protein sequence. "
            "Please check the input sequences."
            "Currently the sequences are: "
            + ", ".join(seq.translate(seq.cdsseq) for seq in self.seqall)
        )
        self.cdsseqall = [seq.cdsseq for seq in self.seqall]
        self.seq_description = exec_options.seq_description
        self.outputdir = exec_options.output
        self.scoringfuncs = scoring_funcs
        self.scoreopts = scoring_options
        self.execopts = exec_options
        self.n_processes = exec_options.processes
        self.quiet = exec_options.quiet
        self.print_top_mutants = exec_options.print_top_mutants
        self.host = exec_options.host
        self.port = exec_options.port
        self.user = exec_options.user
        self.passwd = exec_options.passwd
        self.queue = exec_options.queue

        self.initialize()

    def initialize(self) -> None:
        self.checkpoint_path = os.path.join(self.outputdir, "checkpoints.tsv")
        self.log_file = open(os.path.join(self.outputdir, "log.txt"), "w")
        self.species = self.execopts.species
        self.rand = np.random.RandomState(self.execopts.seed)

        self.mutantall = [
            MutantGenerator(
                seq.rawseq,
                self.rand,
                self.execopts.codon_table,
                self.execopts.protein,
                seq.cdsseq,  # type: ignore[arg-type]
                seq.utr5,
                seq.utr3,
                self.execopts.preserve_stop,
                self.execopts.boost_loop_mutations,
            )
            for seq in self.seqall
        ]
        self.mutantgen = self.mutantall[0]

        if self.execopts.lineardesign_lambda is not None:
            log.info("==> Initializing sequence with LinearDesign...")
            self.mutantgen.lineardesign_initial_codons(
                self.execopts.lineardesign_lambda,
                self.execopts.lineardesign_dir,
                self.execopts.lineardesign_omit_start,
                self.execopts.lineardesign_penalty,
                self.quiet,
            )
        elif self.execopts.random_initialization or self.execopts.protein:
            for mutantgen in self.mutantall:
                mutantgen.randomize_initial_codons()

        if len(self.seqall) == 1 or self.execopts.lineardesign_lambda:
            self.population = [
                [self.mutantgen.utr5]
                + self.mutantgen.initial_codons
                + [self.mutantgen.utr3]
            ]
        else:
            self.population = [
                [seq.utr5] + mutantgen.initial_codons + [seq.utr3]
                for seq, mutantgen in zip(self.seqall, self.mutantall)
            ]

        self.population_foldings: list = [None for _ in self.population]
        self.population_sources: list = [None for _ in self.population]
        parent_no_length = int(np.log10(self.execopts.n_survivors)) + 1
        self.format_parent_no = lambda n, length=parent_no_length: (
            format(n, f"{length}d").replace(" ", "-") if n is not None else "-" * length
        )

        self.mutation_rate = self.execopts.initial_mutation_rate
        self.full_scan_interval = self.execopts.full_scan_interval
        self.in_final_full_scan = False

        self.length_aa = len(self.population[0]) - 2
        self.length_cds = len(self.cdsseqall[0])

        if self.execopts.conservative_start is not None:
            cstart_iter, cstart_width = self.execopts.conservative_start.split(":")
            self.alternative_mutation_fields = [
                (
                    0,
                    int(cstart_iter),
                    self.mutantgen.prepare_alternative_choices(0, int(cstart_width)),
                )
            ]
        else:
            self.alternative_mutation_fields = []
        self.seqeval = SequenceEvaluator(
            self.scoringfuncs,
            self.scoreopts,
            self.execopts,
            self.mutantgen,
            self.species,
            self.length_cds,
            self.quiet,
            self.host,
            self.port,
            self.user,
            self.passwd,
            self.queue,
        )
        self.penalty_metric_flags = self.seqeval.penalty_metric_flags  # XXX
        self.initial_sequence_evaluation = self.seqeval.prepare_evaluation_data(
            "".join(self.population[0])
        )

        self.best_scores: list = []
        self.elapsed_times: list = []
        self.checkpoint_file = open(self.checkpoint_path, "w")
        self.checkpoint_header_written = False

        self.metainfo = {
            "mutation_space": self.mutantgen.compute_mutational_space(),
        }

    def show_configuration(self) -> None:
        spec = self.metainfo["mutation_space"]

        log.info(f"VaxPress Codon Optimizer for mRNA Vaccines {__version__}")
        log.info(hbar_double)
        log.info(f" * Name: {self.seq_description}")
        log.info(f" * CDS length: {self.length_cds} nt")
        log.info(f' * Possible single mutations: {spec["singles"]}')
        log.info(f' * Possible sequences: {spec["total"]}')
        log.info(f' * Command line: {" ".join(sys.argv)}')
        log.info("")

    def next_generation(self, iter_no0: int) -> None:
        # mutate population renamed to next_generation in order to reflect the fact that
        # the population undergoes sexual reproduction
        if (
            self.full_scan_interval > 0
            and (iter_no0 + 1) % self.full_scan_interval == 0
        ):
            return self.prepare_full_scan(iter_no0)

        self.expected_total_mutations = self.mutantgen.compute_expected_mutations(
            self.mutation_rate
        )

        if self.expected_total_mutations < self.stop_threshold:
            if not self.in_final_full_scan:
                log.info(hbar)
                log.info("==> Trying final full scan")
                self.in_final_full_scan = True
                return self.prepare_full_scan(iter_no0)

            log.warning("==> Stopping: expected mutation reaches the minimum")
            raise StopIteration

        self.in_final_full_scan = False
        log.info(hbar)
        log.info(
            f"Iteration {iter_no0+1}/{self.execopts.n_iterations}  -- "
            f"  mut_rate: {self.mutation_rate:.5f} -- "
            f"E(muts): {self.expected_total_mutations:.1f}"
        )

        nextgeneration = self.population[:]
        sources: list[int] = list(range(len(nextgeneration)))

        choices = None
        for begin, end, altchoices in self.alternative_mutation_fields:
            if begin <= iter_no0 < end:
                choices = altchoices
                break

        assert len(self.population) == len(self.population_foldings)

        # population of the next generation
        n_nextgeneration = max(0, self.execopts.n_population - len(self.population))

        for parent_no, _ in zip(
            cycle(list(range(0, len(self.population)))), range(0, n_nextgeneration)
        ):

            parent, parent_folding = (
                self.population[parent_no],
                self.population_foldings[parent_no],
            )
            child = self.mutantgen.generate_mutant(
                parent, self.mutation_rate, choices, parent_folding  # type: ignore[arg-type]
            )
            nextgeneration.append(child)
            sources.append(parent_no)

        self.population[:] = nextgeneration
        self.population_sources[:] = sources
        self.flatten_seqs = ["".join(strand) for strand in self.population]

    def prepare_full_scan(self, iter_no0: int) -> None:
        log.info(hbar)
        log.info(
            f"Iteration {iter_no0+1}/{self.execopts.n_iterations}  -- " "FULL SCAN"
        )

        nextgeneration = self.population[:]
        nextgen_sources = list(range(len(nextgeneration)))

        parents = self.population[:10] if len(self.population) > 10 else self.population
        parents_foldings = (
            self.population_foldings[:10]
            if len(self.population) > 10
            else self.population_foldings
        )

        traverse = self.mutantgen.traverse_all_single_mutations
        for i, (seedseq, seedfold) in enumerate(zip(parents, parents_foldings)):
            for child in traverse(seedseq, seedfold):
                nextgeneration.append(child)
                nextgen_sources.append(i)

        self.population[:] = nextgeneration
        self.population_sources[:] = nextgen_sources
        self.flatten_seqs = ["".join(p) for p in self.population]

    def run(self) -> Iterator[Dict[str, Any]]:
        self.show_configuration()

        timelogs = [time.time()]
        n_survivors = self.execopts.n_survivors
        last_winddown = 0
        error_code = 0

        with futures.ProcessPoolExecutor(max_workers=self.n_processes) as executor:

            if self.execopts.n_iterations == 0:
                # Only the initial sequences are evaluated
                self.flatten_seqs = ["".join(seq) for seq in self.population]
                total_scores, scores, metrics, foldings = self.seqeval.find_fitness(
                    self.flatten_seqs,
                    self.execopts.is_modified,
                    executor,
                    self.execopts.lineardesign_penalty,
                    self.execopts.lineardesign_penalty_weight,
                )
                if total_scores is None:
                    error_code = 1
                else:
                    self.write_checkpoint(
                        0, [0], total_scores, scores, metrics, foldings
                    )
                    timelogs.append(time.time())

            for i in range(self.execopts.n_iterations):
                iter_no = i + 1
                n_parents = len(self.population)

                try:
                    self.next_generation(i)
                except StopIteration:
                    break

                total_scores, scores, metrics, foldings = self.seqeval.find_fitness(
                    self.flatten_seqs,
                    self.execopts.is_modified,
                    executor,
                    self.execopts.lineardesign_penalty,
                    self.execopts.lineardesign_penalty_weight,
                )
                if total_scores is None:
                    # Termination due to errors from one or more scoring functions
                    error_code = 1
                    break

                if len(self.population) > self.execopts.n_population:
                    # Pick the best mutants in each parent to keep diversity
                    idx_sorted = self.prioritized_sort_by_parents(total_scores)
                else:
                    idx_sorted = np.argsort(total_scores)[::-1]
                survivor_indices = idx_sorted[:n_survivors]
                survivors = [self.population[i] for i in survivor_indices]
                survivor_foldings = [foldings[i] for i in survivor_indices]
                self.best_scores.append(total_scores[idx_sorted[0]])

                # Write the evaluation result of the initial sequence in
                # the first iteration
                if i == 0:
                    self.write_checkpoint(
                        0, [0], total_scores, scores, metrics, foldings
                    )

                self.print_eval_results(total_scores, metrics, idx_sorted, n_parents)
                self.write_checkpoint(
                    iter_no, survivor_indices, total_scores, scores, metrics, foldings
                )

                self.population[:] = survivors
                self.population_foldings[:] = survivor_foldings

                log.info(
                    " # Last best scores: "
                    + " ".join(f"{s:.3f}" for s in self.best_scores[-5:])
                )

                if (
                    len(self.best_scores) >= self.execopts.winddown_trigger
                    and iter_no - last_winddown > self.execopts.winddown_trigger
                ):
                    if (
                        self.best_scores[-1]
                        <= self.best_scores[-self.execopts.winddown_trigger]
                    ):
                        self.mutation_rate *= self.execopts.winddown_rate
                        log.info(
                            "==> Winddown triggered: mutation rate "
                            f"= {self.mutation_rate:.5f}"
                        )
                        last_winddown = iter_no

                timelogs.append(time.time())

                self.print_time_estimation(timelogs[-2], timelogs[-1], iter_no)

                log.info("")

                yield {"iter_no": iter_no, "error": error_code, "time": timelogs}

        yield {"iter_no": -1, "error": error_code, "time": timelogs}

    def prioritized_sort_by_parents(self, total_scores):
        bestindices = []

        for src in np.unique(self.population_sources):
            matches = np.where(self.population_sources == src)[0]
            bestidx = max(matches, key=total_scores.__getitem__)
            bestindices.append(bestidx)

        bestindices.sort(key=total_scores.__getitem__, reverse=True)

        # No need to put the rest back because the number of sources equals to
        # the number of survivors.
        return bestindices

    def print_eval_results(self, total_scores, metrics, ind_sorted, n_parents) -> None:
        print_top = min(self.print_top_mutants, len(self.population))
        rowstoshow = ind_sorted[:print_top]
        if len(rowstoshow) < 1:
            return

        metrics_to_show = sorted(
            k
            for k in metrics[rowstoshow[0]].keys()
            if k not in self.penalty_metric_flags
        )
        header = ["flags", "score"] + metrics_to_show
        tabdata = []
        for rank, i in enumerate(rowstoshow):
            flags = [
                (
                    self.format_parent_no(None)
                    if i < n_parents or self.population_sources[i] is None
                    else self.format_parent_no(self.population_sources[i] + 1)
                ),  # parent
                "S " if rank < self.execopts.n_survivors else "- ",
            ]  # is survivor
            for name, flag in self.penalty_metric_flags.items():
                flags.append(flag if metrics[i][name] != 0 in metrics[i] else "-")

            f_total = total_scores[i]
            f_metrics = [metrics[i][name] for name in header[2:]]
            tabdata.append(["".join(flags), f_total] + f_metrics)

        header_short = [h[: self.table_header_length] for h in header]
        log.info(
            tabulate(tabdata, header_short, tablefmt="simple", floatfmt=".2f") + "\n"
        )

    def print_time_estimation(
        self, iteration_start: float, iteration_end: float, iter_no: int
    ) -> None:
        elapsed = iteration_end - iteration_start
        self.elapsed_times.append(elapsed)

        mean_elapsed = np.mean(self.elapsed_times[-5:])
        remaining = (self.execopts.n_iterations - iter_no) * mean_elapsed

        expected_end = time.asctime(time.localtime(time.time() + remaining))
        elapsed_str = time.strftime("%H:%M:%S", time.gmtime(elapsed))
        if elapsed_str.startswith("00:"):
            elapsed_str = elapsed_str[3:]

        log.info(f" # {elapsed_str}s/it  --  Expected finish: {expected_end}")

    def write_checkpoint(
        self, iter_no, survivors, total_scores, scores, metrics, foldings
    ) -> None:
        ind = survivors[0]

        fields = [
            ("iter_no", iter_no),
            ("mutation_rate", self.mutation_rate),
            ("fitness", total_scores[ind]),
        ]
        fields.extend(
            [("metric:" + name, value) for name, value in sorted(metrics[ind].items())]
        )
        fields.extend(
            [("score:" + name, value) for name, value in sorted(scores[ind].items())]
        )
        fields.append(("seq", self.flatten_seqs[ind]))
        fields.append(("structure", foldings[ind]["folding"]))

        if not self.checkpoint_header_written:
            header = [f[0] for f in fields]
            print(*header, sep="\t", file=self.checkpoint_file)
            self.checkpoint_header_written = True

        print(*[f[1] for f in fields], sep="\t", file=self.checkpoint_file)
        self.checkpoint_file.flush()

    def save_results(self):
        # Save the best sequence
        self.bestseq = "".join(self.population[0])
        fastapath = os.path.join(self.outputdir, "best-sequence.fasta")
        with open(fastapath, "w") as f:
            print(f">(best variant)_{self.seq_description}", file=f)
            print(*wrap(self.bestseq, width=self.fasta_line_width), sep="\n", file=f)
            for i in range(1, len(self.population)):
                print(f">(mutant {i})_{self.seq_description}", file=f)
                print(
                    *wrap("".join(self.population[i]), width=self.fasta_line_width),
                    sep="\n",
                    file=f,
                )

        # Save the parameters used for the optimization
        paramspath = os.path.join(self.outputdir, "parameters.json")
        self.save_optimization_parameters(paramspath)

        # Prepare the evaluation results of the best sequence
        return {
            "initial": self.initial_sequence_evaluation,
            "optimized": self.seqeval.prepare_evaluation_data(self.bestseq),
        }

    def save_optimization_parameters(self, path: str) -> None:
        optdata = dump_to_preset(self.scoreopts, self.execopts)
        open(path, "w").write(optdata)
