import ViennaRNA as RNA  # type: ignore

M1PSI_JSON = """\
{
    "modified_base": {
        "name": "N1-Methylpseudouridine (m1Ψ)",
        "sources": [
            {
                "authors": "Mauger, David M., B. Joseph Cabral, Vladimir Presnyak, Stephen V. Su, David W. Reid, Brooke Goodman, Kristian Link et al. ",
                "title": "mRNA structure regulates protein expression through changes in functional half-life",
                "journal": "Proceedings of the National Academy of Sciences 116, no. 48: 24075-24083. ",
                "year": 2019,
                "doi": "10.1073/pnas.1908052116"
            }
        ],
        "unmodified": "U",
        "pairing_partners": [
            "A",
            "G"
        ],
        "one_letter_code": "1",
        "fallback": "U",
        "stacking_energies": {
            "G1CA": -2.43,
            "1CAG": -2.67,
            "C1GA": -1.83,
            "1GAC": -2.26,
            "A11A": -1.13,
            "11AA": -1.18,
            "1AA1": -1.86
        },
        "terminal_energies": {
            "1A": 0.31,
            "A1": 0.31
        }
    }
}
"""

M1PSI_PARAMS = RNA.sc_mod_read_from_json(M1PSI_JSON)

if __name__ == "__main__":
    new_seq = "AGGACAUUUGCUUCUGACACAACUGUGUUCACUAGCAACCUCAAACAGACACCAUGGCCGUUUACCCAUACGAUGUUCCUGACUAUGCGGGCUAUCCCUAUGACGUCCCGGACUAUGCAGGCUCCUAUCCAUAUGACGUUCCAGAUUACGCUGGAUCUGGCGUCUUCACACUCGAAGAUUUCGUUGGGGACUGGCGACAGACAGCCGGCUACAACCUGGACCAAGUCCUUGAACAGGGAGGUGUGUCCAGUUUGUUUCAGAAUCUCGGGGUGUCCGUAACUCCGAUCCAAAGGAUUGUCCUGAGCGGUGAAAAUGGGCUGAAGAUCGACAUCCAUGUCAUCAUCCCGUAUGAAGGUCUGAGCGGCGACCAAAUGGGCCAGAUCGAAAAAAUUUUUAAGGUGGUGUACCCUGUGGAUGAUCAUCACUUUAAGGUGAUCCUGCACUAUGGCACACUGGUAAUCGACGGGGUUACGCCGAACAUGAUCGACUAUUUCGGACGGCCGUAUGAAGGCAUCGCCGUGUUCGACGGCAAAAAGAUCACUGUAACAGGGACCCUGUGGAACGGCAACAAAAUUAUCGACGAGCGCCUGAUCAACCCCGACGGCUCCCUGCUGUUCCGAGUAACCAUCAACGGAGUGACCGGCUGGCGGCUGUGCGAACGCAUUCUGGCGUAAGCUCGCUUUCUUGCUGUCCAAUUUCUAUUAAAGGUUCCUUUGUUCCCUAAGUCCAACUACUAAACUGGGGGAUAUUAUGAAGGGCCUUGAGCAUCUGGAUUCUGCCUAAUAAAAAACAUUUAUUUUCAUUGC"
    fc = RNA.fold_compound(new_seq)
    fc.sc_mod(
        params=M1PSI_PARAMS,
        modification_sites=[
            i for i in range(1, len(new_seq) + 1) if new_seq[i - 1] == "U"
        ],
    )
    struct, efe = fc.pf()
    mea_struct = fc.MEA()
    print(f"m1Ψ energy: {efe}, MEA structure: {mea_struct}")
    print(f"else structure is {struct}")
