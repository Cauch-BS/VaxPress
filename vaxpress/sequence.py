# added at 2024-06-20
# defined the Sequence object which the mutant generate object inherits

from typing import Union

import pandas as pd  # type: ignore[import-untyped]
from Bio.Data import CodonTable
from vaxpress.log import log

PROTEIN_ALPHABETS = "ACDEFGHIKLMNPQRSTVWY" + "*"
RNA_ALPHABETS = "ACGU"

STOP = "*"


class Sequence:
    def __init__(
        self,
        rawseq: Union[str, list[str]],
        codon_table: str = "standard",
        is_protein: bool = False,
        cdsseq: str = "",
        utr5: str = "",
        utr3: str = "",
        preserve_stop: bool = False,
    ):
        self.initialize_codon_table(codon_table)
        self.codons: list = []

        self.utr5: str = utr5
        self.utr3: str = utr3
        self.rawseq = rawseq

        if is_protein:
            protein_seq = rawseq or cdsseq
            invalid_letters = set(protein_seq) - set(PROTEIN_ALPHABETS)
            if invalid_letters:
                raise ValueError(
                    "Invalid protein sequence: " f'{" ".join(invalid_letters)}'
                )
            if isinstance(protein_seq, str) and protein_seq[-1] != "*":
                protein_seq += "*"
            elif (
                isinstance(protein_seq, list)
                and isinstance(protein_seq[0], str)
                and protein_seq[-1] != "*"
            ):
                protein_seq.append("*")
            self.cdsseq: str = self.backtranslate(protein_seq)

        else:
            if isinstance(rawseq, str):
                self.rawseq = rawseq.upper().replace("T", "U")
            elif isinstance(rawseq, list) and isinstance(rawseq[0], str):
                self.rawseq = "".join(rawseq).upper().replace("T", "U")

            invalid_letters = set(self.rawseq) - set(RNA_ALPHABETS)
            if invalid_letters:
                raise ValueError(
                    "Invalid RNA sequence: " f'{" ".join(invalid_letters)}'
                )

        if not cdsseq and not is_protein:
            self._preserve_stop = preserve_stop
            self.utr5, self.cdsseq, self.codons, self.utr3 = self.truncate(
                self.rawseq, self._preserve_stop
            )
        elif cdsseq and isinstance(cdsseq, str) and not is_protein:
            self.cdsseq = cdsseq
            self.utr5 = utr5
            self.utr3 = utr3
            self.rawseq = utr5 + cdsseq + utr3
            self.codons = [
                self.cdsseq[i : i + 3] for i in range(0, len(self.cdsseq), 3)
            ]
        elif is_protein:
            self.rawseq = protein_seq
            self.utr5 = utr5
            self.utr3 = utr3
            self.codons = [
                self.cdsseq[i : i + 3] for i in range(0, len(self.cdsseq), 3)
            ]
        else:
            raise ValueError("Invalid CDS sequence!")

        if len(self.cdsseq) % 3 != 0:
            raise ValueError("Invalid CDS sequence length!")

    @staticmethod
    def truncate(
        rawseq: str | list[str], preserve_stop: bool = False
    ) -> tuple[str, str, list[str], str]:
        if isinstance(rawseq, list):
            rawseq = "".join(rawseq)
        # Find the first AUG codon with a Kozak sequence
        start = rawseq.find("AUG")
        if start == -1:
            raise ValueError("CDS start not found in sequence")
        new_start, i = start, 0
        while rawseq[new_start - 3] not in {"A", "G"}:
            new_start = rawseq.find("AUG", new_start + 1)
            if new_start == -1:
                new_start = start
                break
            i += 1
        start = new_start
        # make a reading frame from the start codon
        reading_frame = [rawseq[i : i + 3] for i in range(start, len(rawseq), 3)]
        min_stop_index = len(reading_frame)
        stop_codons = ["UAA", "UAG", "UGA"]
        for stop in stop_codons:
            try:
                stop_index = reading_frame.index(stop)
                if min_stop_index > stop_index:
                    min_stop_index = stop_index
            except ValueError:
                continue
        if not preserve_stop:
            reading_frame_length = (
                min_stop_index + 1
                if min_stop_index < len(reading_frame)
                else len(reading_frame)
            )
        elif preserve_stop:
            reading_frame_length = min_stop_index
        reading_frame = reading_frame[:reading_frame_length]
        utr5 = rawseq[:start]
        utr3 = rawseq[start + 3 * (reading_frame_length) :]
        return utr5, "".join(reading_frame), reading_frame, utr3

    def initialize_codon_table(self, codon_table: str) -> None:
        table_var_name = f"{codon_table}_rna_table"
        if not hasattr(CodonTable, table_var_name):
            raise ValueError(f"Invalid codon table name: {codon_table}")

        self.codon_table = getattr(CodonTable, table_var_name)

        codon_frame = pd.DataFrame(
            list(self.codon_table.forward_table.items())
            + [[stopcodon, STOP] for stopcodon in self.codon_table.stop_codons],
            columns=["codon", "aa"],
        )

        self.synonymous_codons, self.aa2codons, self.codon2aa = {}, {}, {}

        for aa, codons in codon_frame.groupby("aa"):
            codons = set(codons["codon"])
            self.aa2codons[aa] = codons

            for codon in codons:
                self.synonymous_codons[codon] = sorted(codons - set([codon]))
                self.codon2aa[codon] = aa

    def backtranslate(self, proteinseq: Union[str, list]) -> str:
        return "".join(next(iter(self.aa2codons[aa])) for aa in proteinseq)

    def translate(self, rnaseq: str) -> str:
        return "".join(
            self.codon2aa[rnaseq[i : i + 3]] for i in range(0, len(rnaseq), 3)
        )

    def get_polyA(self) -> int:
        """Returns the (reversed) index of first polyA tail in the 3' UTR.
        If no polyA tail is found, return 1"""
        length = len(self.utr3)
        idx = self.utr3.find("AAAAAAAA")
        result = -idx if idx == -1 else length - idx
        return result


if __name__ == "__main__":
    seq = "GCCACCATGATTTAAAAAAAAAAAAAAAAAAAAAA"
    cdsseq = Sequence(seq)
    print(cdsseq.codons)
    print(cdsseq.utr5)
    print(cdsseq.utr3)
    print(seq[: -cdsseq.get_polyA()])
    protein_seq = "MTHRL*"
    protein = Sequence(protein_seq, is_protein=True)
    print(protein.rawseq)
