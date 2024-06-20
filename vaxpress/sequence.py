#added at 2024-06-20
#defined the Sequence object which the mutant generate object inherits

from Bio.Data import CodonTable
import pandas as pd
from typing import Union

STOP = "*"

class Sequence:
    def __init__(self, cdsseq: Union[str, list], codon_table: str = 'standard', is_protein: bool = False):
        self.initialize_codon_table(codon_table)

        if is_protein:
            self.cdsseq = self.backtranslate(cdsseq)
        else:
            if type(cdsseq) == str:
                self.cdsseq = cdsseq.upper()
            elif type(cdsseq) == list and type(cdsseq[0]) == str:
                self.cdsseq = ''.join(cdsseq)


        if len(self.cdsseq) % 3 != 0:
            raise ValueError("Invalid CDS sequence length!")

        self.aaseq = self.translate(self.cdsseq)
        self.codons = [self.cdsseq[i:i+3] for i in range(0, len(self.cdsseq), 3)]

    def initialize_codon_table(self, codon_table: str) -> None:
        table_var_name = f'{codon_table}_rna_table'
        if not hasattr(CodonTable, table_var_name):
            raise ValueError(f'Invalid codon table name: {codon_table}')

        self.codon_table = getattr(CodonTable, table_var_name)

        codon_frame = pd.DataFrame(
            list(self.codon_table.forward_table.items()) +
            [[stopcodon, STOP] for stopcodon in self.codon_table.stop_codons],
            columns = ['codon', 'aa']
        )

        self.synonymous_codons, self.aa2codons, self.codon2aa = {}, {}, {}

        for aa, codons in codon_frame.groupby('aa'):
            codons = set(codons['codon'])
            self.aa2codons[aa] = codons

            for codon in codons:
                self.synonymous_codons[codon] = sorted(codons - set([codon]))
                self.codon2aa[codon] = aa

    def backtranslate(self, proteinseq: str) -> str:
        return ''.join(next(iter(self.aa2codons[aa])) for aa in proteinseq)

    def translate(self, rnaseq: str) -> str:
        return ''.join(self.codon2aa[rnaseq[i:i+3]] for i in range(0, len(rnaseq), 3))