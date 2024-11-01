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

import numpy as np

from vaxpress.scoring import ScoringFunction


class PairingProbFitness(ScoringFunction):

    name = "aup"
    description = "average unpaired probability"
    priority = 101
    uses_folding = False
    uses_basepairing_prob = True

    arguments = [
        (
            "weight",
            dict(
                metavar="WEIGHT",
                type=float,
                default=-5.0,
                help="scoring weight for AUP (default: -5.0)",
            ),
        ),
        (
            "u-weight",
            dict(
                metavar="U_WEIGHT",
                type=float,
                default=3.0,
                help="weight for U (default: 3.0)",
            ),
        ),
        (
            "a-weight",
            dict(
                metavar="A_WEIGHT",
                type=float,
                default=1.5,
                help="weight for A (default: 1.5)",
            ),
        ),
    ]

    def __init__(self, weight, u_weight, a_weight, _length_cds):
        self.weight = weight
        self.u_weight = u_weight
        self.a_weight = a_weight
        self.length = _length_cds
        self.trans_table = bytes.maketrans(b"ACGU", b"\x00\x01\x02\x03")

    def score(self, seqs, foldings):
        weighted_aups = []
        weights = np.array([self.a_weight, 1, 1, self.u_weight])  # ACGU
        for seq, folding in zip(seqs, foldings):
            pi_array = folding["pi_array"]
            seqindex = np.frombuffer(
                seq.encode().translate(self.trans_table), dtype=np.uint8
            )
            # map weights to index
            wi = np.choose(seqindex, weights)
            xi = 1.0 - pi_array
            weighted_aup = np.dot(xi, wi) / np.sum(wi)
            weighted_aup_val = weighted_aup.item()
            weighted_aups.append(weighted_aup_val)
        scores = [weighted_aup * self.weight for weighted_aup in weighted_aups]
        return {self.name: scores}, {self.name: weighted_aups}

    def annotate_sequence(self, seq, folding):
        weights = np.array([self.a_weight, 1, 1, self.u_weight])
        pi_array = folding["pi_array"]
        seqindex = np.frombuffer(
            seq.encode().translate(self.trans_table), dtype=np.uint8
        )
        wi = np.choose(seqindex, weights)
        xi = 1 - pi_array
        weighted_aup = np.dot(xi, wi) / np.sum(wi)
        weighted_aup_val = weighted_aup.item()
        return {"aup": weighted_aup_val}

    def evaluate_local(self, seq, folding):
        xi_array = 1 - folding["pi_array"]
        baseindex = list(range(len(seq)))
        return {"gc": (baseindex, xi_array)}
