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

from collections import namedtuple
from typing import Iterable, Union

import numpy as np

from . import lineardesign
from .sequence import Sequence

MutationChoice = namedtuple("MutationChoice", ["pos", "altcodon"])


class MutantGenerator(Sequence):

    rawseq: Union[str, list[str]] = []
    initial_codons: list = []

    def __init__(
        self,
        rawseq: Union[str, list],
        random_state: np.random.RandomState,
        codon_table: str = "standard",
        is_protein: bool = False,
        cdsseq: str = "",
        utr5: str = "",
        utr3: str = "",
        preserve_stop: bool = False,
        boost_loop_mutations: str = None,  # type: ignore[assignment]
    ):
        super().__init__(
            rawseq, codon_table, is_protein, cdsseq, utr5, utr3, preserve_stop
        )
        self.rand = np.random if random_state is None else random_state

        if boost_loop_mutations is not None:
            tokens = boost_loop_mutations.split(":")
            (self.boost_loop_mutations_weight, self.boost_loop_mutations_start) = float(
                tokens[0]
            ), int(tokens[1])

        else:
            self.boost_loop_mutations_weight = 0

        self.initial_codons = self.codons[:]
        self.setup_choices()

    def setup_choices(self) -> None:
        choices = []
        initial_codons = []

        for i in range(len(self.cdsseq) // 3):
            codon = self.cdsseq[i * 3 : i * 3 + 3]
            alternatives = len(self.synonymous_codons[codon])
            for alt in range(alternatives):
                choices.append(MutationChoice(i, alt))
            initial_codons.append(codon)

        self.choices = choices
        self.initial_codons = initial_codons

    def prepare_alternative_choices(self, left, right) -> list[MutationChoice]:
        return [choice for choice in self.choices if left <= choice.pos < right]

    def randomize_initial_codons(self) -> None:
        self.initial_codons[:] = [
            self.rand.choice([codon] + self.synonymous_codons[codon])
            for codon in self.initial_codons
        ]

    def lineardesign_initial_codons(
        self, lmd, lddir, omitstart, penalty_region, quiet
    ) -> None:
        prot = self.translate(self.cdsseq)

        prot_om = prot[omitstart:]
        res = lineardesign.run_lineardesign(lddir, prot_om, penalty_region, lmd, quiet)
        assert len(res["seq"]) == len(res["str"]) == len(prot_om) * 3

        rseq = res["seq"]
        self.initial_codons[omitstart:] = [
            rseq[i * 3 : i * 3 + 3] for i in range(len(prot_om))
        ]

    def calc_probabilities(
        self, choices: list[MutationChoice], folding: dict
    ) -> np.ndarray:
        minimum_position = self.boost_loop_mutations_start
        loop_positions = set(
            [i // 3 for i, code in enumerate(folding["folding"]) if code == "."]
        )

        weightmap = [1, self.boost_loop_mutations_weight]
        probs = [
            weightmap[mut.pos in loop_positions and mut.pos >= minimum_position]
            for mut in choices
        ]
        return np.array(probs) / sum(probs)

    def generate_mutant(
        self,
        codons: list[str],
        mutation_rate: float,
        choices: list[MutationChoice] = None,  # type: ignore[assignment]
        folding: dict = None,  # type: ignore[assignment]
    ) -> list[str]:
        child = codons[:]
        if choices is None:
            choices = self.choices

        # Draw number of mutations from binomial distribution
        n_mutations = self.rand.binomial(len(choices), mutation_rate)
        n_mutations = max(1, min(n_mutations, len(choices)))

        prob_dist = (
            self.calc_probabilities(choices, folding)
            if self.boost_loop_mutations_weight != 0 and folding is not None
            else None
        )

        # Select mutations
        mutation_choices = self.rand.choice(
            len(choices), n_mutations, replace=False, p=prob_dist
        )

        # Apply mutations
        for i in mutation_choices:
            mut = choices[i]
            child[mut.pos + 1] = self.synonymous_codons[child[mut.pos + 1]][
                mut.altcodon
            ]

        return child

    def traverse_all_single_mutations(
        self,
        parent: list[str],
        fold: dict[str, Iterable],
        choices: list[MutationChoice] = None,  # type: ignore[assignment]
    ):
        if choices is None:
            choices = self.choices

        loop_positions = set(
            [i // 3 for i, code in enumerate(fold["folding"]) if code == "."]
        )

        for choice in self.choices:
            if choice.pos not in loop_positions:
                continue

            child = parent[:]
            child[choice.pos + 1] = self.synonymous_codons[child[choice.pos + 1]][
                choice.altcodon
            ]
            yield child

    def compute_expected_mutations(self, mutation_rate: float) -> float:
        return len(self.choices) * mutation_rate

    def compute_mutational_space(self) -> dict:
        log10_totalcases = np.sum(np.log10(list((map(len, self.choices)))))
        totalcases_mantissa = 10 ** (log10_totalcases % 1)
        totalcases = f"{totalcases_mantissa:.2f} x 10^{int(log10_totalcases)}"

        return {
            "singles": len(self.choices),
            "total": totalcases,
        }
