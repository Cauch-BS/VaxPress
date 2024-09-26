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

from vaxpress.scoring import ScoringFunction


class UnpairedUridineFitness(ScoringFunction):

    name = "unpaired_u"
    description = "average unpaired uridine count"
    priority = 102
    uses_folding = False
    uses_basepairing_prob = True

    arguments = [
        (
            "weight",
            dict(
                metavar="WEIGHT",
                type=float,
                default=-3.0,
                help="weight for unpaired ucount (default: -3.0)",
            ),
        ),
    ]

    def __init__(self, weight, _length_cds):
        self.weight = weight
        self.length = _length_cds

    def score(self, seqs, pairingprobs):
        U_sups = []
        for seq, pairingprob in zip(seqs, pairingprobs):
            U_idx = [i for i, base in enumerate(seq) if base == "U"]
            pi_cooarray = pairingprob["pi_array"]
            pi_array = pi_cooarray.sum(axis=0)
            # sum only the Pi of the U index
            total_unpairedu_probs = sum(
                1 - pi_array[i] for i in U_idx
            )  # to be minimized
            U_sups.append(total_unpairedu_probs)
        scores = [u_sup * self.weight for u_sup in U_sups]

        return {self.name: scores}, {self.name: U_sups}

    def annotate_sequence(self, seq, pairingprob):
        U_idx = [i for i, base in enumerate(seq) if base == "U"]
        pi_cooarray = pairingprob["pi_array"]
        pi_array = pi_cooarray.sum(axis=0)
        U_sup = sum(1 - pi_array[i] for i in U_idx)
        return {self.name: U_sup}
