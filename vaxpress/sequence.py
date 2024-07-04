#added at 2024-06-20
#defined the Sequence object which the mutant generate object inherits

from Bio.Data import CodonTable
import pandas as pd
from typing import Union

STOP = "*"

class Sequence:
    def __init__(self, cdsseq: Union[str, list], codon_table: str = 'standard', 
                 is_protein: bool = False, is_cds = False):
        self.initialize_codon_table(codon_table)
        self.codons = []
        self.utr5 = ''
        self.utr3 = ''

        if is_protein:
            self.cdsseq = self.backtranslate(cdsseq)

        else:
            if type(cdsseq) == str:
                self.cdsseq = cdsseq.upper().replace('T', 'U')
            elif type(cdsseq) == list and type(cdsseq[0]) == str:
                self.cdsseq = ''.join(cdsseq)

        if not is_cds:
            self.utr5, self.cdsseq, self.codons, self.utr3 = self.truncate(self.cdsseq)

        if len(self.cdsseq) % 3 != 0:
            raise ValueError("Invalid CDS sequence length!")
        
        if not self.codons:
            self.codons = [self.cdsseq[i:i+3] for i in range(0, len(self.cdsseq), 3)]
    
    @staticmethod
    def truncate(rawseq: str) -> str:
        start = rawseq.find('AUG')
        reading_frame = [rawseq[i:i+3] for i in range(start, len(rawseq), 3)]
        min_stop_index = len(reading_frame)
        stop_codons = ['UAA', 'UAG', 'UGA']
        for stop in stop_codons:
            try:
                stop_index = reading_frame.index(stop)
                if min_stop_index > stop_index:
                    min_stop_index = stop_index
            except ValueError:
                continue
        reading_frame = reading_frame[:(min_stop_index + 1)]
        utr5 = rawseq[:start]
        utr3 = rawseq[start + 3*(min_stop_index + 1):]
        return utr5, ''.join(reading_frame), reading_frame, utr3

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
    
    def get_polyA(self) -> int:
        '''Returns the (reversed) index of first polyA tail in the 3' UTR.
        If no polyA tail is found, return 1'''
        length = len(self.utr3)
        idx = self.utr3.find('AAAAAAAA')
        result = - idx if idx == -1 else length - idx
        return result
    
if __name__ == "__main__":
    seq = 'CCCCCATGATTTAAAAAAAAAAAAAAAAAAAAAA'
    cdsseq = Sequence(seq)
    print(cdsseq.codons)
    print(cdsseq.utr5)
    print(cdsseq.utr3)
    print(seq[:- cdsseq.get_polyA()])