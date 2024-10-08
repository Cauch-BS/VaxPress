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
            pi_cooarray = folding["pi_array"]
            pi_array = pi_cooarray.sum(axis=0)
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
        pi_cooarray = folding["pi_array"]
        pi_array = pi_cooarray.sum(axis=0)
        xi_array = 1 - pi_array
        U_sup = np.dot(xi_array, U_idx)
        U_sup = U_sup.item()
        return {self.name: U_sup}
