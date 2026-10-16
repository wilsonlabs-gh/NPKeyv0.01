# NPKey

Annotated genome -> antiSMASH BGCs -> BiG-SLiCE GCF assignment -> reference compounds -> SQLite/JSON/CSV.

| File | Responsibility |
|---|---|
| `main.py` | Call the five stages. |
| `npkey/ingest.py` | Stage 1: preserve and validate local GenBank or a versioned NCBI assembly. |
| `npkey/antismash.py` | Stage 2: submit once, resume polling, download or import antiSMASH results. |
| `npkey/gcf_assignment.py` | Stage 3: query the frozen reference and apply the assignment cutoff. |
| `npkey/compound_context.py` | Stage 4: link reference BGCs to documented compounds and resolve InChIKeys. |
| `npkey/save_results.py` | Stage 5: save evidence JSON, the SQLite ledger and a CSV summary. |
| `npkey/gcf_reference.py` | Build/cache MIBiG-derived reference families once. |
| `npkey/antismash_results.py` | Verify input sequences and region files; normalize coordinates. |
| `npkey/io_utils.py` | Shared checksums, atomic JSON, read retries, SQL reads and command logging. |

See [QCinput.md](QCinput.md) for genome requirements. `tests/` holds the offline
tests and their fixtures. Keep local inputs/reference databases in `data/` and
outputs in `runs/`; both directories and the virtual environments are ignored by Git.

## Install

The pipeline uses Python 3.12+. BiG-SLiCE runs separately on Linux with Python 3.11;
its older Biopython dependency must not replace the pipeline's dependencies.
Use Linux or WSL for the complete workflow. Create both environments in this folder:

```sh
python3.12 -m venv .venv
.venv/bin/pip install -r requirements.txt
python3.11 -m venv .bigslice
.bigslice/bin/pip install -r requirements-bigslice.txt
.bigslice/bin/download_bigslice_hmmdb
```

## Prepare the reference once

Download the [official antiSMASH-processed MIBiG 4.0 GenBank archive](https://dl.secondarymetabolites.org/mibig/mibig_antismash_4.0_gbk_as8b1.tar.bz2)
to `data/`. Also download the matching [MIBiG 4.0 JSON metadata archive](https://dl.secondarymetabolites.org/mibig/mibig_json_4.0.tar.gz)
for stage 4. Family preparation requires the processed GenBank archive.
Run these examples with the pipeline Python environment:

```python
from npkey.gcf_reference import prepare_reference

reference = prepare_reference(
    "data/mibig_antismash_4.0_gbk_as8b1.tar.bz2", "data/reference",
    ".bigslice/bin/bigslice-models", executable=".bigslice/bin/bigslice",
    threshold=0.4,
)
```

The cache records the release/checksum, engine/HMM versions, clustering settings,
family members and database checksum. Repeating preparation reuses a completed
reference or resumes its unfinished build. Changed inputs require a new directory.
Family IDs belong to this model; a new build gets a new model identifier.

## Run a genome

Save this as a Python script in the project root and run it with
`.venv/bin/python your_script.py`:

```python
from main import run

result = run("data/genome.gbk", "runs/sample", "data/reference",
             mibig_metadata="data/mibig_json_4.0.tar.gz",
             executable=".bigslice/bin/bigslice")
print(result)
# For an existing antiSMASH 8 result, also pass saved_results="data/antismash_output".
```

Stage 2 preserves annotations, uses relaxed detection and disables KnownClusterBlast.
It saves the job ID before polling; repeating a completed run reuses its files.
Uncertain submissions require recovering the existing job ID. Completed stage-2
JSON survives downstream failures. Imported results must match the original genome
and include region GenBank files. Coordinates are zero-based and end-exclusive;
origin-crossing regions retain their ordered segments.

Stage 3 returns ranked GCF candidates and distances. The nearest candidate is
accepted only within `cutoff` (default: reference clustering threshold). This is a
configurable distance rule, not a probability or experimentally validated identity.
Unassigned regions retain their nearest candidates. Zero detected features is
reported as `no_features`; zero BGCs is valid. No assignment against MIBiG does not
establish novelty.

Queries reuse completed local results; interrupted local queries can be rerun.
Each query copies the database and centroid cache, leaving the reference frozen.
`provenance.gcf` identifies the model and its MIBiG family members.

Stage 4 reads the matching MIBiG metadata release for accepted family members.
Explicit PubChem CIDs resolve to full InChIKeys; identical keys share a compound
node. Missing, failed or conflicting identifiers retain unresolved nodes with null
InChIKeys. Cached PubChem responses record retrieval time, URL and checksum;
the ledger also records the metadata archive checksum and lookup errors.

SQLite and evidence JSON retain reference BGC accessions, compound nodes, documented
BGC-compound links and region-to-reference links labelled `same_gcf`. These are
reference associations; they do not establish which compounds a query BGC produces.
CSV summarizes GCF assignments.
`result` contains the run ID and JSON/CSV/SQLite paths. Repeated identical saves
reuse rows rather than duplicating them. Raw antiSMASH files and BiG-SLiCE query
databases/logs remain available under the run directory.

## Verify

```sh
.venv/bin/python -m unittest discover -s tests -v
```

Tests cover offline handoff, reference caching, assignment boundaries, failed jobs,
replay, compound identity, provenance and persistence. BiG-SLiCE commands are mocked;
a [published-data fixture](tests/fixtures/bigslice_v2/README.md) also checks real
version-2 scores in the upstream query schema. A live HMM/clustering run remains untested.
Sources: [antiSMASH API](https://api.secondarymetabolites.org/),
[BiG-SLiCE query mode](https://github.com/medema-group/bigslice#querying-antismash-bgcs),
[MIBiG releases](https://dl.secondarymetabolites.org/mibig/).
