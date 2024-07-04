from . import ScoringFunction
import ctypes 
import os
repeat_path = os.path.abspath('./repeats.so')
repeats = ctypes.CDLL(repeat_path) 
repeats.returnRepeatsPenalty.argtypes =  [ctypes.c_char_p, ctypes.c_int]
repeats.returnRepeatsPenalty.restype = ctypes.c_float

class RepeatsFitness(ScoringFunction):
    name = 'repeats'
    description = 'Remote Repeats'
    priority = 35

    arguments = [
        ('weight',
         dict(type=float, default=1.0, metavar='WEIGHT',
              help='scoring weight for remote repeats (default: 1.0)')),
        ('min-length',
         dict(type=int, default=8, metavar='LENGTH',
              help='minimum length of repeats to be considered as a remote '
                   'repeat (default: 8)')),
    ]

    penalty_metric_flags = {'remote repeat': 'rr'}

    def __init__(self, weight, min_length, _length_cds):
        self.weight = weight 
        self.min_length = min_length

    def score(self, seqs):
        penalties = []
        for seq in seqs:
            seq = seq.encode('utf-8')
            penalty = repeats.returnRepeatsPenalty(seq, self.min_length)
            penalties.append(penalty)

        remote_repeat_score = [penalty * self.weight for penalty in penalties]

        metrics = {'repeat': penalties}
        scores = {'repeat': remote_repeat_score}

        return scores, metrics
