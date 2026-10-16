"""Offline tests using synthetic GenBank records and mocked NCBI responses."""

import hashlib
from io import BytesIO, StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock, patch
import requests
from zipfile import ZipFile

from Bio import BiopythonParserWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from ingest import ingest_genome


ACCESSION = "GCF_000203835.1"
MEMBER = f"ncbi_dataset/data/{ACCESSION}/genomic.gbff"


def record(identifier="chromosome", sequence="ATGAAATAA", cds=True):
    """Invented nine-base fixture, not a biological reference sequence."""
    result = SeqRecord(Seq(sequence), id=identifier, name=identifier)
    result.annotations["molecule_type"] = "DNA"
    if cds:
        result.features = [SeqFeature(SimpleLocation(0, 9, strand=1), type="CDS")]
    return result


def genbank(*records):
    handle = StringIO()
    SeqIO.write(records, handle, "genbank")
    return handle.getvalue().encode()


def package(contents, member=MEMBER):
    handle = BytesIO()
    with ZipFile(handle, "w") as archive:
        archive.writestr(member, contents)
        archive.writestr("../../unwanted.txt", "Must never be extracted")
    return handle.getvalue()


class IngestTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.output = self.root / "run"
        self.source = self.root / "input.gbff"

    def local(self, contents):
        self.source.write_bytes(contents)
        return ingest_genome(self.source, self.output)

    def test_preserves_all_records_and_bytes_and_returns_stage_two_input(self):
        contents = genbank(record(), record("plasmid"), record("small_contig", cds=False))
        result = self.local(contents)
        saved = Path(result["genbank_path"])
        self.assertEqual(set(result), {"genbank_path", "accession", "checksum"})
        self.assertTrue(saved.is_absolute())
        self.assertEqual(saved.read_bytes(), contents)
        self.assertEqual(self.source.read_bytes(), contents)
        self.assertEqual(result["checksum"], hashlib.sha256(contents).hexdigest())
        self.assertIsNone(result["accession"])
        with saved.open() as handle:
            self.assertEqual([r.id for r in SeqIO.parse(handle, "genbank")],
                             ["chromosome", "plasmid", "small_contig"])

    def test_repeated_ingest_and_reingesting_saved_file_are_identical(self):
        first = self.local(genbank(record()))
        self.assertEqual(first, ingest_genome(self.source, self.output))
        self.assertEqual(first, ingest_genome(first["genbank_path"], self.output))
        self.assertEqual(len(list((self.output / "input").glob("*.gbk"))), 1)

    def test_formatting_warning_is_visible_but_does_not_reject_or_rewrite_input(self):
        contents = genbank(record()).replace(b"ORIGIN\n", b"ORIGIN\n\n")
        with self.assertWarnsRegex(BiopythonParserWarning, "Blank line in sequence data"):
            result = self.local(contents)
        self.assertEqual(Path(result["genbank_path"]).read_bytes(), contents)
        self.assertEqual(result["checksum"], hashlib.sha256(contents).hexdigest())

    def test_rejects_invalid_inputs_without_leaving_artifacts(self):
        invalid = [
            (b"", "No GenBank"),
            (b">sequence\nATGAAATAA\n", "No GenBank"),
            (genbank(record(), record()), "Duplicate sequence"),
            (genbank(record(cds=False)), "No CDS annotations"),
            (genbank(record(sequence="ATGZZZTAA")), "IUPAC DNA"),
            (genbank(record(sequence="", cds=False)), "IUPAC DNA"),
            (genbank(record()).replace(b"9 bp", b"99 bp"), "Invalid annotated"),
            (genbank(record()).replace(b"//\n", b""), "Invalid annotated"),
        ]
        for contents, message in invalid:
            with self.subTest(message=message, contents=contents[:40]):
                with self.assertRaisesRegex(ValueError, message):
                    self.local(contents)
                self.assertEqual(list((self.output / "input").iterdir()), [])

    def test_rejects_out_of_bounds_cds(self):
        invalid = record()
        invalid.features[0].location = SimpleLocation(0, 12)
        with self.assertRaisesRegex(ValueError, "out-of-bounds"):
            self.local(genbank(invalid))

    def test_keeps_origin_crossing_cds_unchanged(self):
        circular = record()
        circular.annotations["topology"] = "circular"
        circular.features[0].location = SimpleLocation(6, 9) + SimpleLocation(0, 3)
        contents = genbank(circular)
        result = self.local(contents)
        self.assertEqual(Path(result["genbank_path"]).read_bytes(), contents)

    def test_rejects_undefined_sequence(self):
        undefined = record()
        undefined.seq = Seq(None, length=9)
        with self.assertRaisesRegex(ValueError, "Invalid annotated"):
            self.local(genbank(undefined))

    def test_rejects_protein_record(self):
        protein = record()
        protein.annotations["molecule_type"] = "protein"
        with self.assertRaisesRegex(ValueError, "DNA GenBank"):
            self.local(genbank(protein))

    def test_missing_file_and_unversioned_accession_do_not_call_network(self):
        with patch("io_utils.requests.get") as network:
            with self.assertRaises(FileNotFoundError):
                ingest_genome(self.source, self.output)
            with self.assertRaisesRegex(ValueError, "versioned assembly"):
                ingest_genome("GCF_000203835", self.output)
            network.assert_not_called()

    def test_download_preserves_exact_assembly_and_does_not_extract_other_files(self):
        contents = genbank(record(), record("plasmid"))
        with patch("io_utils.requests.get", return_value=Mock(content=package(contents))) as network:
            result = ingest_genome(ACCESSION, self.output)
        self.assertEqual(result["accession"], ACCESSION)
        self.assertEqual(Path(result["genbank_path"]).read_bytes(), contents)
        url = network.call_args.args[0]
        self.assertIn(f"/{ACCESSION}/download?include_annotation_type=GENOME_GBFF", url)
        self.assertEqual(list(self.root.rglob("unwanted.txt")), [])

    def test_does_not_substitute_another_assembly_or_missing_annotation(self):
        for member in (MEMBER.replace(".1/", ".2/"), MEMBER.replace("gbff", "fna")):
            with self.subTest(member=member):
                with patch("io_utils.requests.get", return_value=Mock(content=package(genbank(record()), member))):
                    with self.assertRaisesRegex(ValueError, "no unique annotated"):
                        ingest_genome(ACCESSION, self.output)
                self.assertEqual(list((self.output / "input").iterdir()), [])

    def test_invalid_zip_is_explicit(self):
        with patch("io_utils.requests.get", return_value=Mock(content=b"not a zip")):
            with self.assertRaisesRegex(ValueError, "invalid or damaged ZIP"):
                ingest_genome(ACCESSION, self.output)

    def test_retries_transient_reads(self):
        responses = [requests.ConnectionError("temporary failure"), Mock(content=package(genbank(record())))]
        with patch("io_utils.requests.get", side_effect=responses) as network, patch("io_utils.time.sleep") as sleep:
            ingest_genome(ACCESSION, self.output)
        self.assertEqual(network.call_count, 2)
        sleep.assert_called_once_with(1)

    def test_permanent_http_error_is_not_retried(self):
        response = requests.Response()
        response.status_code = 404
        error = requests.HTTPError(response=response)
        with patch("io_utils.requests.get", side_effect=error) as network:
            with self.assertRaises(requests.RequestException):
                ingest_genome(ACCESSION, self.output)
        self.assertEqual(network.call_count, 1)

    def test_exhausted_retries_leave_no_artifact(self):
        with patch("io_utils.requests.get", side_effect=requests.ConnectionError("offline")) as network, patch("io_utils.time.sleep"):
            with self.assertRaises(requests.RequestException):
                ingest_genome(ACCESSION, self.output)
        self.assertEqual(network.call_count, 3)
        self.assertEqual(list((self.output / "input").iterdir()), [])


if __name__ == "__main__":
    unittest.main()
