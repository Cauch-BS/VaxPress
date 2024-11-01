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

from ..sequence import Sequence
from . import ScoringFunction


class OverallRepeatsFitness(ScoringFunction):

    name = "orepeats"
    description = "Overall Repeats"
    priority = 60

    arguments = [
        (
            "weight",
            dict(
                type=float,
                default=10.0,
                metavar="WEIGHT",
                help="scoring weight for overall repeats (default: 1.0)",
            ),
        ),
        (
            "min-repeat-length",
            dict(
                type=int,
                default=8,
                metavar="N",
                help="minimum length of a repeats to be considered as a "
                "repeat (default: 2)",
            ),
        ),
    ]

    penalty_metric_flags = {"orepeat": "or"}

    def __init__(self, weight, min_repeats, min_length, _length_cds):
        self.weight = weight
        self.min_repeats = min_repeats
        self.min_length = min_length
        self.reverse = str.maketrans("ACGT", "TGCA")

    def reverse_complement(self, seq):
        return seq.translate(self.reverse)[::-1]

    def scan_repeats(self, seq, length):
        repeatpos = set()
        for i in range(len(seq) - length + 1):
            if i in repeatpos:
                continue
            fwdseq = seq[i : i + length]
            fwdseq_rc = self.reverse_complement(fwdseq)
            for j in range(i + 1, len(seq)):
                if j in repeatpos:
                    continue
                revseq = seq[j : j + length]
                if revseq == fwdseq:
                    repeatpos.add(i)
                    repeatpos.add(j)
                elif seq[j : j + length] == fwdseq_rc:
                    repeatpos.add(i)
                    repeatpos.add(j)
        repmarks = np.zeros(len(seq) - length + 1, dtype=int)
        for i in repeatpos:
            repmarks[i : i + length] = 1
        return repmarks.mean()

    def score(self, seqs):
        rep_percent = []
        for seq in seqs:
            seq = seq[: -Sequence(seq).get_polyA()]
            repeats = self.scan_repeats(seq, self.min_length)
            rep_percent.append(repeats)

        repeat_score = [per * self.weight for per in rep_percent]

        metrics = {"orepeat": rep_percent}
        scores = {"orepeat": repeat_score}

        return scores, metrics
