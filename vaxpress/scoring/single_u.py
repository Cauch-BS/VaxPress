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
import numpy as np


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
        self.trans_table = bytes.maketrans(b"ACGU", b"\x00\x01\x02\x03")

    def score(self, seqs, foldings):
        U_sups = []
        for seq, folding in zip(seqs, foldings):
            seqindex = np.frombuffer(
                seq.encode().translate(self.trans_table), dtype=np.uint8
            )
            u_idx = np.where(seqindex == 3, 1, 0)
            pi_array = folding["pi_array"]
            xi_array = 1 - pi_array
            # sum only the Pi of the U index
            total_unpaired_u_probs = np.dot(xi_array, u_idx)
            unpaired_u_prob_val = total_unpaired_u_probs.item()
            U_sups.append(unpaired_u_prob_val)
        scores = [u_sup * self.weight for u_sup in U_sups]

        return {self.name: scores}, {self.name: U_sups}

    def annotate_sequence(self, seq, folding):
        seqindex = np.frombuffer(
            seq.encode().translate(self.trans_table), dtype=np.uint8
        )
        U_idx = np.where(seqindex == 3, 1, 0)
        pi_array = folding["pi_array"]
        xi_array = 1 - pi_array
        U_sup = np.dot(xi_array, U_idx)
        U_sup = U_sup.item()
        return {self.name: U_sup}
