"""Controller: ingest, mine, assign families, add compounds, persist."""

from npkey.ingest import ingest_genome
from npkey.antismash import run_antismash
from npkey.gcf_assignment import assign_gcfs
from npkey.compound_context import add_compounds
from npkey.save_results import save_results


def run(source, output_dir, reference_dir, saved_results=None, *, mibig_metadata,
        cutoff=None, executable="bigslice"):
    genome = ingest_genome(source, output_dir)
    evidence = run_antismash(genome, output_dir, saved_results=saved_results)
    evidence = assign_gcfs(evidence, reference_dir, output_dir, cutoff=cutoff, executable=executable)
    evidence = add_compounds(evidence, mibig_metadata, output_dir)
    return save_results(evidence, output_dir)
