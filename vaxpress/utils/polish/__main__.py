import click
from Bio.SeqIO import parse, to_dict

from .polish_helper import (
    load_sequence_structure,
    find_alternatives,
    make_structure,
    write_fasta,
)
from .idt_complexity import IDTComplexity


@click.group()
def polish():
    """CLI tool for sequence polishing and complexity analysis."""
    pass


@polish.command()
@click.option("-i", "--input", type=click.File("r"), required=True)
def find_complexity(input):
    def read_as_dna(input):
        seqs = to_dict(parse(input, "fasta"))
        return {
            str(seq_id): str(seq.seq).upper().replace("U", "T")
            for seq_id, seq in seqs.items()
        }

    seqs = read_as_dna(input)
    scores = IDTComplexity().score(seqs.keys(), seqs.values())
    click.echo("Complexity scores: \n")
    for name, scores in scores.items():
        click.echo(f"{name} : {scores['Total Score']}")
        click.echo(f"Violated rules: {scores['Violated']}")


@polish.command()
@click.option(
    "--cds-start",
    type=int,
    required=False,
    default=-1,
    help="Start position of CDS (in 1-based index)",
)
@click.option(
    "--cds-end",
    type=int,
    required=False,
    default=-1,
    help="End position of CDS (in 1-based index)",
)
@click.option(
    "-o", "--output", type=click.File("w"), default="codon_polish_results.fasta"
)
@click.option("-i", "--input", type=click.File("r"), required=True)
@click.option(
    "--roi-start",
    type=int,
    required=True,
    help="Start position of ROI (in 1-based index)",
)
@click.option(
    "--roi-end",
    type=int,
    required=True,
    help="End position of ROI (in 1-based index)",
)
def identify_alternatives(output, input, cds_start, cds_end, roi_start, roi_end):
    try:
        seq, structure = load_sequence_structure(input)
    except ValueError:
        seq, structure = make_structure(input)
    if cds_start == -1:
        click.echo("CDS start not provided. Attempting to find CDS start automatically")
        cds_start = seq.find("AUG")
        if cds_start == -1:
            raise ValueError("CDS start not found in sequence")
        i = 0
        new_cds_start = cds_start
        while seq[new_cds_start - 3] not in {"A", "G"} and new_cds_start != -1:
            click.echo(
                f"Invalid Kozak sequence. Attempting to find new CDS start: Attempt {i}"
            )
            new_cds_start = seq.find("AUG", new_cds_start + 1)  # Find next AUG
            i += 1
        cds_start = new_cds_start

        reading_frame = [seq[i : i + 3] for i in range(cds_start, len(seq), 3)]
        min_stop_index = len(reading_frame)
        stop_codons = ["UAA", "UAG", "UGA"]
        for stop in stop_codons:
            try:
                stop_index = reading_frame.index(stop)
                if min_stop_index > stop_index:
                    min_stop_index = stop_index
            except ValueError:
                continue
        valid_cds = reading_frame[: (min_stop_index + 1)]
        cds_end = cds_start + 3 * len(valid_cds)

    roi_start -= 1

    if cds_start < 0 or roi_start < 0:
        raise ValueError("Invalid start position. ROI start must be in 1-based index")

    if (cds_end - cds_start) % 3 != 0:
        raise ValueError("Invalid CDS length")

    if roi_start < cds_start or roi_end > cds_end:
        raise ValueError("Region of interest is outside of CDS")

    combinations = find_alternatives(
        seq, structure, cds_start, cds_end, roi_start, roi_end
    )
    write_fasta(output, seq, structure, cds_start, cds_end, combinations)
    click.echo(f"Results written to {output.name}")
