#
# VaxPress
#
# Copyright 2023-2024 Seoul National University
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

from ...data.codon_usage_data import codon_usage
from ...data import codon_table
from itertools import groupby
import numpy as np
import pandas as pd
import RNA

def setup_lookup_tables():
    synonymous_codons = {}
    aa2codons = {}

    sorted_codons = sorted(codon2aa.items(), key=lambda x: x[::-1])
    for aa, codons in groupby(sorted_codons, key=lambda x: x[1]):
        codons = [codon for codon, aa in codons]
        aa2codons[aa] = codons

        for c in codons:
            synonymous_codons[c] = sorted(set(codons) - {c})
    
    return synonymous_codons, aa2codons

codon_usage = codon_usage['Homo sapiens']
codon2aa = codon_table.standard_table
synonymous_codons, aa2codons = setup_lookup_tables()

def make_structure(input):
    input.seek(0)
    seq = ''
    for line in input:
        if line.startswith('>'):
            continue
        seq_line = line.strip().upper().replace('T','U')
        seq += seq_line
    folded = RNA.fold_compound(seq)
    (ss, _) = folded.mfe()
    return seq, ss

def load_sequence_structure(input):
    input.seek(0)
    seq = []
    structure = []
    for line in input:
        if line.startswith('>'):
            continue

        seq_entry = line.strip().upper().replace('T', 'U')
        try:
            structure_entry = next(input).strip()
        except StopIteration:
            raise ValueError('Structure unavailable for some lines')

        if set(structure_entry) - set('.()'):
            raise ValueError('Invalid structure: {}'.format(structure_entry))

        seq.append(seq_entry)
        structure.append(structure_entry)

    return ''.join(seq), ''.join(structure)

def find_alternatives(seq, structure, cds_start, cds_end, roi_start, roi_end):
    print(f'ROI Sequence:  {seq[roi_start:roi_end]}')
    print(f'ROI Structure: {structure[roi_start:roi_end]}')
    col_seqbegin = len('ROI Structure: ')

    iter_roi_structure = enumerate(structure[roi_start:roi_end], roi_start)
    unpaired_positions = [i for i, s in iter_roi_structure if s == '.']

    combinations = []
    for i in unpaired_positions:
        frame = (i - cds_start) % 3
        codon_start = i - frame
        codon = seq[codon_start:codon_start+3]
        origcodon_usage = codon_usage[codon]
        alternatives = synonymous_codons[codon]

        for altcodon in alternatives:
            mutations = [j for j, (orig, alt) in enumerate(zip(codon, altcodon))
                         if orig != alt]
            if frame not in mutations:
                continue

            altcodon_usage = codon_usage[altcodon]
            combinations.append((codon_start, codon, altcodon,
                                 np.log2(altcodon_usage / origcodon_usage)))

    for codon_start, codon, altcodon, deltacai in combinations:
        relposition = codon_start - roi_start
        print(str(codon_start + 1).rjust(col_seqbegin + relposition - 1), end=' ')

        altcodon_emphasized = ''.join(
            b.lower() if a == b
            else (f'\x1b[1;36m{b}\x1b[0m' if
                  (f + codon_start) in unpaired_positions else f'\x1b[1;31m{b}\x1b[0m')
        for f, (a, b) in enumerate(zip(codon, altcodon)))
        print(altcodon_emphasized, end=' ')

        print(f'  {deltacai:.2f}')
    
    return combinations

def add_frame_annotation(cdsseq):
    cdsseq = list(cdsseq.lower())
    for i in range(0, len(cdsseq), 3):
        cdsseq[i] = cdsseq[i].upper()
    return ''.join(cdsseq)

def write_fasta(output, seq, structure, cds_start, cds_end, combinations):
    annotated_seq = (seq[:cds_start] + add_frame_annotation(seq[cds_start:cds_end].lower()) +
                     seq[cds_end:])

    for codon_start, codon, altcodon, deltacai in combinations:
        output.write(f'>{codon_start + 1}:{codon}->{altcodon} deltaCAI:{deltacai:.2f}\n')

        altseq = annotated_seq[:codon_start] + altcodon + annotated_seq[codon_start+3:]

        output.write(altseq + '\n')
        output.write(structure + '\n')