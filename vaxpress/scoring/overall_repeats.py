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

from . import ScoringFunction
from ..sequence import Sequence
from collections import defaultdict


def find_kmer_repeats_percentage(rna_seq, min_length=8):
    kmer_count = defaultdict(int)
    total_repeats = len(rna_seq) - min_length + 1

    for i in range(total_repeats):
        forward_kmer = rna_seq[i : i + min_length]
        kmer_count[forward_kmer] += 1
        reverse_kmer = forward_kmer[::-1].translate(str.maketrans("AUCG", "UAGC"))
        kmer_count[reverse_kmer] += 1

    total_repeated_kmers = sum([v for v in kmer_count.values() if v > 1])
    percentage_repeats = total_repeated_kmers / total_repeats
    return percentage_repeats


class OverallRepeatsFitness(ScoringFunction):

    name = "o_repeats"
    description = "Overall Repeats"
    priority = 61

    arguments = [
        (
            "weight",
            dict(
                type=float,
                default=5.0,
                metavar="WEIGHT",
                help="scoring weight for overall repeats (default: 5.0)",
            ),
        ),
        (
            "min-length",
            dict(
                type=int,
                default=2,
                metavar="LENGTH",
                help="minimum length of k-mer repeats to be considered"
                "repeat (default: 8)",
            ),
        ),
    ]

    penalty_metric_flags = {"o_repeat": "or"}

    def __init__(self, weight, min_length, _length_cds):
        self.weight = weight
        self.min_length = min_length

    def score(self, seqs):
        replengths = []
        for seq in seqs:
            seq = seq[: -Sequence(seq).get_polyA()]
            repeat_metrics = find_kmer_repeats_percentage(seq, self.min_length)
            replengths.append(repeat_metrics)

        repeat_score = [metric * self.weight for metric in replengths]

        metrics = {"o_repeat": replengths}
        scores = {"o_repeat": repeat_score}

        return scores, metrics
