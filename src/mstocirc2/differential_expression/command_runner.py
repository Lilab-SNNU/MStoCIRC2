
                "  mstocirc2 dea -pm peptide_matrix.tsv -cr circ_translation_evidence.tsv -de design.txt -st generic -o out_dea",
    p.add_argument("-cr", "--circrna-reference", metavar="<FILE>", required=True, help="Path to the circRNA reference table, for example 'circ_translation_evidence.tsv'.")