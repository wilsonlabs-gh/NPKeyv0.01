"""Stage 1: annotated GenBank or versioned NCBI accession -> preserved genome."""

from io import BytesIO
from pathlib import Path
import re
import shutil
from tempfile import TemporaryDirectory
import warnings
from zipfile import BadZipFile, ZipFile

from Bio import BiopythonParserWarning, SeqIO
from .io_utils import checksum, get


def ingest_genome(source, output_dir):
    source_path = Path(source).expanduser()
    accession = None
    if not source_path.is_file():
        if re.fullmatch(r"GC[AF]_[0-9]{9}\.[1-9][0-9]*", str(source)):
            accession = str(source)
        elif str(source).startswith(("GCF_", "GCA_")):
            raise ValueError("Use a versioned assembly accession, e.g. GCF_000203835.1.")
        else:
            raise FileNotFoundError(source_path)
    input_dir = Path(output_dir).expanduser().resolve() / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".ingest-", dir=input_dir) as temporary:
        staged = Path(temporary) / "genome.gbk"
        if accession:
            _download_genbank(accession, staged)
        else:
            shutil.copyfile(source_path, staged)
        _validate_genbank(staged)
        digest = checksum(staged)
        destination = input_dir / f"{digest}.gbk"
        staged.replace(destination)
    return {"genbank_path": str(destination), "accession": accession, "checksum": digest}


def _download_genbank(accession, destination):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"
    content = get(url + "?include_annotation_type=GENOME_GBFF").content
    member = f"ncbi_dataset/data/{accession}/genomic.gbff"
    try:
        with ZipFile(BytesIO(content)) as archive:
            if archive.namelist().count(member) != 1:
                raise ValueError(f"NCBI package has no unique annotated GenBank file for {accession}.")
            with archive.open(member) as source, destination.open("wb") as handle:
                shutil.copyfileobj(source, handle)
    except BadZipFile as exc:
        raise ValueError(f"NCBI returned an invalid or damaged ZIP for {accession}.") from exc


def _validate_genbank(path):
    identifiers = set()
    cds_count = 0
    try:
        with path.open(encoding="utf-8") as handle, warnings.catch_warnings():
            warnings.filterwarnings("error", "Expected sequence length", BiopythonParserWarning)
            warnings.filterwarnings("error", "Premature end of file", BiopythonParserWarning)
            for record in SeqIO.parse(handle, "genbank"):
                if not record.id or record.id == "<unknown id>" or any(c.isspace() for c in record.id):
                    raise ValueError("Each sequence needs an identifier without whitespace.")
                if record.id in identifiers:
                    raise ValueError(f"Duplicate sequence identifier: {record.id}")
                identifiers.add(record.id)
                sequence = str(record.seq).upper()
                if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
                    raise ValueError(f"{record.id}: expected a non-empty, defined IUPAC DNA sequence.")
                if "DNA" not in record.annotations.get("molecule_type", "").upper():
                    raise ValueError(f"{record.id}: expected a DNA GenBank record.")
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    if feature.location is None:
                        raise ValueError(f"{record.id}: CDS has no readable location.")
                    for part in feature.location.parts:
                        if part.ref is not None:
                            raise ValueError(f"{record.id}: CDS refers to an external sequence.")
                        start, end = int(part.start), int(part.end)
                        if start < 0 or end > len(sequence) or start >= end:
                            raise ValueError(f"{record.id}: CDS location is empty or out-of-bounds.")
                    cds_count += 1
    except (ValueError, TypeError, BiopythonParserWarning) as exc:
        raise ValueError(f"Invalid annotated GenBank input: {exc}") from exc
    if not identifiers:
        raise ValueError("No GenBank sequence records found; supply annotated GenBank, not FASTA.")
    if not cds_count:
        raise ValueError("No CDS annotations found; supply an annotated GenBank genome.")
