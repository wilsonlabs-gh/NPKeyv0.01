"""Invented sequences in the antiSMASH 8 JSON/GenBank format; no biological claims."""

from io import BytesIO, StringIO
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock, patch
from zipfile import ZipFile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
import requests

from npkey.antismash import SETTINGS, run_antismash
from npkey.antismash_results import parse_results
from npkey.ingest import ingest_genome
from npkey.io_utils import get


def result_data():
    return {"version": "8.0.4", "schema": 4, "records": [{
        "id": "TEST.1", "seq": {"data": "ATG" * 40},
        "areas": [{"start": 3, "end": 90, "products": ["NRPS"]}],
        "features": [{"type": "region", "qualifiers": {
            "region_number": ["1"], "contig_edge": ["False"], "rules": ["test-rule"]}}],
        "modules": {"antismash.detection.hmm_detection": {
            "strictness": "relaxed", "enabled_types": ["NRPS"],
            "rule_results": {"multipliers": {"cutoff": 1.0}}}},
    }]}


def region_files(data):
    files = {}
    for record in data["records"]:
        for number, area in enumerate(record["areas"], 1):
            sequence = record["seq"]["data"]
            start, end = area["start"], area["end"]
            sequence = sequence[start:end] if start < end else sequence[start:] + sequence[:end]
            bgc = SeqRecord(Seq(sequence), id=record["id"], name=record["id"])
            bgc.annotations = {"molecule_type": "DNA", "structured_comment": {
                "antiSMASH-Data": {"Version": "8.0.4"}}}
            loc = SimpleLocation(0, len(sequence), strand=1)
            bgc.features = [SeqFeature(loc, type="region", qualifiers={
                "region_number": [str(number)], "contig_edge": ["False"], "product": ["NRPS"]}),
                SeqFeature(loc, type="protocluster"),
                SeqFeature(loc, type="CDS", qualifiers={"translation": ["M" * (len(sequence) // 3)]})]
            handle = StringIO()
            SeqIO.write(bgc, handle, "genbank")
            files[f"{record['id']}.region{number:03d}.gbk"] = handle.getvalue()
    return files


def response(*, data=None, text="", content=b""):
    result = Mock(content=content, text=text)
    result.json.return_value = data
    return result


def archive_bytes(data, unsafe=False):
    handle = BytesIO()
    with ZipFile(handle, "w") as archive:
        archive.writestr("genome.json", json.dumps(data))
        for name, contents in region_files(data).items():
            archive.writestr(name, contents)
        if unsafe:
            archive.writestr("../escape.txt", "unsafe")
    return handle.getvalue()


class StageTwoTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.output = self.root / "run"
        source = self.root / "input.gbk"
        record = SeqRecord(Seq("ATG" * 40), id="TEST.1", name="TEST")
        record.annotations["molecule_type"] = "DNA"
        record.features = [SeqFeature(SimpleLocation(0, 120), type="CDS")]
        SeqIO.write(record, source, "genbank")
        self.genome = ingest_genome(source, self.output)
        self.saved = self.root / "saved"
        self.saved.mkdir()
        self.data = result_data()
        self.job_dir = self.output / "antismash" / self.genome["checksum"]

    def write(self):
        (self.saved / "genome.json").write_text(json.dumps(self.data))
        for name, contents in region_files(self.data).items():
            (self.saved / name).write_text(contents)

    def parse(self):
        self.write()
        return parse_results(self.saved, self.genome)

    def reads(self, unsafe=False):
        return [response(data={"state": "done", "result_url": "/upload/job/index.html"}),
                response(text='<a href="genome.zip">Download</a>'),
                response(content=archive_bytes(self.data, unsafe))]

    def checkpoint(self):
        self.job_dir.mkdir(parents=True, exist_ok=True)
        (self.job_dir / "job.json").write_text(json.dumps({"job_id": "job-1", "parameters": SETTINGS}))

    def test_handoff_preserves_coordinates_ids_and_real_bgc_paths(self):
        evidence = self.parse()
        region = evidence["regions"][0]
        self.assertEqual((region["start"], region["end"]), (3, 90))
        self.assertEqual(region["segments"], [[3, 90]])
        self.assertFalse(region["contig_edge"])
        self.assertTrue(Path(region["bgc_path"]).is_file())
        self.assertEqual(len(region["bgc_checksum"]), 64)
        self.assertEqual(set(evidence), {"regions", "gcf_matches", "issues", "provenance"})
        self.assertEqual(evidence, self.parse())

    def test_region_id_tracks_detection_settings(self):
        first = self.parse()["regions"][0]["region_id"]
        self.data["records"][0]["modules"]["antismash.detection.hmm_detection"]["strictness"] = "strict"
        self.assertNotEqual(first, self.parse()["regions"][0]["region_id"])

    def test_origin_crossing_region_retains_segments_and_sequence(self):
        self.data["records"][0]["areas"][0].update(start=100, end=20)
        region = self.parse()["regions"][0]
        self.assertEqual(region["segments"], [[100, 120], [0, 20]])
        self.assertTrue(region["crosses_origin"])

    def test_zero_regions_and_obsolete_module_are_ignored(self):
        self.data["records"][0]["modules"]["antismash.modules.clusterblast"] = "ignored legacy content"
        self.assertEqual(self.parse()["issues"], [])
        self.data["records"][0]["areas"] = []
        self.assertEqual(self.parse()["regions"], [])

    def test_rejects_genome_mismatch_and_missing_records(self):
        self.data["records"][0]["seq"]["data"] = "A" * 120
        with self.assertRaisesRegex(ValueError, "differs from the input"):
            self.parse()
        self.data["records"] = []
        with self.assertRaisesRegex(ValueError, "every input sequence"):
            self.parse()

    def test_renamed_record_uses_original_sequence_id_and_output_filename(self):
        self.data["records"][0].update(id="renamed", original_id="TEST.1")
        region = self.parse()["regions"][0]
        self.assertEqual(region["sequence_id"], "TEST.1")
        self.assertEqual(Path(region["bgc_path"]).name, "renamed.region001.gbk")

    def test_rejects_modified_input(self):
        with Path(self.genome["genbank_path"]).open("a") as handle:
            handle.write("\n")
        with self.assertRaisesRegex(ValueError, "checksum"):
            self.parse()

    def test_rejects_missing_and_wrong_region_files(self):
        self.write()
        bgc = self.saved / "TEST.1.region001.gbk"
        bgc.unlink()
        with self.assertRaisesRegex(ValueError, "region GenBank"):
            parse_results(self.saved, self.genome)
        self.write()
        bgc.write_text(bgc.read_text().replace("atgatgatga", "aaaaaaaaaa"))
        with self.assertRaisesRegex(ValueError, "BGC sequence"):
            parse_results(self.saved, self.genome)

    def test_rejects_unsupported_and_malformed_json(self):
        self.data["schema"] = 999
        with self.assertRaisesRegex(ValueError, "Unsupported antiSMASH"):
            self.parse()
        self.data["schema"] = 4
        del self.data["records"][0]["modules"]
        with self.assertRaisesRegex(ValueError, "Malformed antiSMASH"):
            self.parse()

    def test_import_bypasses_network_and_preserves_stage_two_evidence(self):
        self.write()
        with patch("npkey.antismash.requests.post") as post, patch("npkey.antismash.get") as get_mock:
            evidence = run_antismash(self.genome, self.output, saved_results=self.saved)
        post.assert_not_called()
        get_mock.assert_not_called()
        self.assertEqual(json.loads((self.job_dir / "evidence.json").read_text()), evidence)

    def test_submit_download_and_replay(self):
        with patch("npkey.antismash.requests.post", return_value=response(data={"id": "job-1"})) as post:
            with patch("npkey.antismash.get", side_effect=self.reads()):
                evidence = run_antismash(self.genome, self.output)
        self.assertEqual(post.call_count, 1)
        self.assertEqual(post.call_args.kwargs["data"]["knownclusterblast"], "false")
        with patch("npkey.antismash.requests.post") as post, patch("npkey.antismash.get") as get_mock:
            self.assertEqual(run_antismash(self.genome, self.output), evidence)
        post.assert_not_called()
        get_mock.assert_not_called()

    def test_timeout_and_partial_download_resume_same_job(self):
        self.checkpoint()
        with patch("npkey.antismash.get", return_value=response(data={"state": "running"})):
            with self.assertRaises(TimeoutError):
                run_antismash(self.genome, self.output, timeout=0)
        (self.job_dir / "results").mkdir()
        with patch("npkey.antismash.requests.post") as post, patch("npkey.antismash.get", side_effect=self.reads()):
            self.assertEqual(len(run_antismash(self.genome, self.output)["regions"]), 1)
        post.assert_not_called()

    def test_uncertain_submission_is_not_repeated(self):
        with patch("npkey.antismash.requests.post", side_effect=requests.Timeout()) as post:
            with self.assertRaisesRegex(RuntimeError, "not confirmed"):
                run_antismash(self.genome, self.output)
            with self.assertRaisesRegex(RuntimeError, "uncertain"):
                run_antismash(self.genome, self.output)
        self.assertEqual(post.call_count, 1)

    def test_failed_job_and_changed_settings_do_not_resubmit(self):
        self.checkpoint()
        with patch("npkey.antismash.requests.post") as post:
            with patch("npkey.antismash.get", return_value=response(data={"state": "failed"})):
                with self.assertRaisesRegex(RuntimeError, "failed"):
                    run_antismash(self.genome, self.output)
            with patch("npkey.antismash.SETTINGS", dict(SETTINGS, genefinder="different")):
                with self.assertRaisesRegex(ValueError, "Settings changed"):
                    run_antismash(self.genome, self.output)
        post.assert_not_called()

    def test_zip_traversal_is_rejected(self):
        self.checkpoint()
        with patch("npkey.antismash.get", side_effect=self.reads(unsafe=True)):
            with self.assertRaisesRegex(ValueError, "outside its directory"):
                run_antismash(self.genome, self.output)
        self.assertFalse((self.job_dir / "escape.txt").exists())

    def test_read_retries_only_transient_failures(self):
        good = response()
        with patch("npkey.io_utils.requests.get", side_effect=[requests.Timeout(), good]) as network:
            with patch("npkey.io_utils.time.sleep"):
                self.assertIs(get("https://example.invalid"), good)
        self.assertEqual(network.call_count, 2)
        failed = requests.Response()
        failed.status_code = 404
        with patch("npkey.io_utils.requests.get", return_value=failed) as network:
            with self.assertRaises(requests.HTTPError):
                get("https://example.invalid")
        self.assertEqual(network.call_count, 1)
