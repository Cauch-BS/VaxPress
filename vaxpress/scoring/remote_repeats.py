from . import ScoringFunction
from ..sequence import Sequence
import repeats 


class RepeatsFitness(ScoringFunction):
    name = 'repeats'
    description = 'Remote Repeats'
    priority = 35

    arguments = [
        ('weight',
         dict(type=float, default=5.0, metavar='WEIGHT',
              help='scoring weight for remote repeats (default: 5.0)')),
        ('min-length',
         dict(type=int, default=8, metavar='LENGTH',
              help='minimum length of repeats to be considered as a remote '
                   'repeat (default: 8)')),
    ]

    penalty_metric_flags = {'remote_repeat': 'rr'}

    def __init__(self, weight, min_length, _length_cds):
        self.weight = weight 
        self.min_length = min_length

    def score(self, seqs):
        penalties = []
        for seq in seqs:
            seq = seq[:- Sequence(seq).get_polyA()]
            print(seq[-10:])
            penalty = repeats.returnRepeatsPenalty(seq, self.min_length)
            penalties.append(penalty)

        remote_repeat_score = [penalty * self.weight for penalty in penalties]

        metrics = {'remote_repeat': penalties}
        scores = {'remote_repeat': remote_repeat_score}

        return scores, metrics
