
        help=(
            "Path to the circRNA ORF mapping table. Optional; if omitted together with "
            "-fi/--circ-file-input, -co/--circ-orf is used as both the truncated and full-length ORF."
        ),
    if not args.circ_file_input and not args.circ_mapping:
        log.warning(
            "Neither -fi/--circ-file-input nor -cm/--circ-mapping was provided; "
            "using -co/--circ-orf as both the BSJ-truncated ORF and the full-length source ORF."
        )