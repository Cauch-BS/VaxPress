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
            str(seq_id): str(seq.seq).replace("U", "T") for seq_id, seq in seqs.items()
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
    required=True,
    default=1,
    help="Start position of CDS (in 1-based index)",
)
@click.option(
    "--cds-end", type=int, required=True, help="End position of CDS (in 1-based index)"
)
@click.option(
    "-o", "--output", type=click.File("w"), default="codon_polish_results.fasta"
)
@click.option("-i", "--input", type=click.File("r"), required=True)
@click.option("--roi-start", type=int, required=True)
@click.option("--roi-end", type=int, required=True)
def identify_alternatives(cds_start, cds_end, output, input, roi_start, roi_end):
    cds_start -= 1
    roi_start -= 1

    if cds_start < 0 or roi_start < 0:
        raise ValueError("Invalid start position")

    if (cds_end - cds_start) % 3 != 0:
        raise ValueError("Invalid CDS length")

    if roi_start < cds_start or roi_end > cds_end:
        raise ValueError("Region of interest is outside of CDS")

    try:
        seq, structure = load_sequence_structure(input)
    except ValueError:
        seq, structure = make_structure(input)

    combinations = find_alternatives(
        seq, structure, cds_start, cds_end, roi_start, roi_end
    )
    write_fasta(output, seq, structure, cds_start, cds_end, combinations)
    click.echo(f"Results written to {output.name}")
