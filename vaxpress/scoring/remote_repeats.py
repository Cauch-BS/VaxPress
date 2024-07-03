
from . import ScoringFunction

class RepeatsFitness(ScoringFunction):

    name = 'repeats'
    description = 'Remote Repeats'
    priority = 25

    arguments = [
        ('weight',
         dict(type=float, default=1.0, metavar='WEIGHT',
              help='scoring weight for remote repeats (default: 1.0)')),
        ('min-repeats',
         dict(type=int, default=2, metavar='N',
              help='minimum number of repeats to be considered as a remote '
                   'repeat (default: 2)')),
        ('min-length',
         dict(type=int, default=10, metavar='LENGTH',
              help='minimum length of repeats to be considered as a remote '
                   'repeat (default: 8)')),
    ]

    penalty_metric_flags = {'remote repeat': 'rr'}

    def __init__(self, weight, min_repeats, min_length, _length_cds):
        self.weight = weight / _length_cds * -1000
        self.min_repeats = min_repeats
        self.min_length = min_length

    def score(self, seqs):
        replengths = []
        for seq in seqs:
            repeats = pytrf.GTRFinder('name', seq, min_repeat=self.min_repeats,
                                    min_length=self.min_length)
            replengths.append(sum(r.length for r in repeats))

        repeat_score = [length * self.weight for length in replengths]

        metrics = {'repeat': replengths}
        scores = {'repeat': repeat_score}

        return scores, metrics
