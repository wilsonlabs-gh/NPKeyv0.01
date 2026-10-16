"""Reference chemistry, identifier resolution and ledger integration."""

from contextlib import closing
from copy import deepcopy
from io import BytesIO
import json
from pathlib import Path
import sqlite3
import tarfile
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock, patch
import requests

from compound_context import add_compounds
from io_utils import checksum
from save_results import save_results

KEY = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
OTHER_KEY = "RYYVLZVUVIJVGH-UHFFFAOYSA-N"


def write_metadata(path, records, release="4.0"):
    with tarfile.open(path, "w:gz") as archive:
        for accession, record in records.items():
            content = json.dumps(record).encode()
            member = tarfile.TarInfo(f"mibig_json_{release}/{accession}.json")
            member.size = len(content)
            archive.addfile(member, BytesIO(content))


class CompoundTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.archive = self.root / "mibig.tar.gz"
        self.records = {
            "BGC0000001": {"accession": "BGC0000001", "version": 5, "quality": "questionable",
                           "status": "active", "compounds": [
                               {"name": "first label", "structure": "CC", "databaseIds": ["pubchem:2244"]},
                               {"name": "unresolved", "databaseIds": []}]},
            "BGC0000002": {"accession": "BGC0000002", "version": 3, "compounds": [
                {"name": "second label", "databaseIds": ["pubchem:2244"]},
                {"name": "unresolved", "databaseIds": []}]}}
        self.evidence = {"regions": [
            {"region_id": "r1", "sequence_id": "seq", "start": 0, "end": 100,
             "gcf_id": 10, "gcf_status": "assigned"},
            {"region_id": "r2", "sequence_id": "seq", "start": 100, "end": 200,
             "gcf_id": None, "gcf_status": "unassigned"}], "gcf_matches": [], "issues": [],
            "provenance": {"gcf": {"mibig_release": "4.0", "model_id": "fixture-model", "cutoff": 0.4,
                "families": {"10": ["BGC0000001", "BGC0000002"], "20": ["BGC9999999"]}}}}
        self.sleep = patch("compound_context.time.sleep")
        self.sleep.start()
        self.addCleanup(self.sleep.stop)
        write_metadata(self.archive, self.records)

    def response(self, cid=2244, key=KEY):
        response = Mock()
        response.json.return_value = {"PropertyTable": {"Properties": [{"CID": cid, "InChIKey": key}]}}
        return response

    def enrich(self):
        return add_compounds(deepcopy(self.evidence), self.archive, self.root)

    def test_explicit_cid_deduplication_preserves_source_labels_and_unresolved_nodes(self):
        with patch("compound_context.get", return_value=self.response()) as get:
            result = self.enrich()
        get.assert_called_once_with("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/InChIKey/JSON")
        self.assertEqual(len(result["compounds"]), 3)
        self.assertEqual(result["compounds"][KEY], KEY)
        self.assertEqual(sum(value is None for value in result["compounds"].values()), 2)
        self.assertEqual(result["reference_bgcs"]["BGC0000001"]["version"], 5)
        self.assertEqual(result["region_bgcs"], [("r1", "BGC0000001"), ("r1", "BGC0000002")])
        for accession, identity, index, compound in result["bgc_compounds"]:
            self.assertEqual(compound, self.records[accession]["compounds"][index])
        provenance = result["provenance"]["compounds"]
        self.assertEqual(provenance["sha256"], checksum(self.archive))
        lookup = provenance["pubchem"]["2244"]
        self.assertEqual(lookup["sha256"], checksum(lookup["snapshot"]))
        self.assertIn("retrieved_at", lookup)

    def test_cached_replay_and_ledger_links_are_repeatable(self):
        with patch("compound_context.get", return_value=self.response()):
            first = self.enrich()
        with patch("compound_context.get", side_effect=AssertionError("Must use cache")):
            self.assertEqual(self.enrich(), first)
        saved = save_results(first, self.root)
        self.assertEqual(saved, save_results(self.enrich(), self.root))
        with closing(sqlite3.connect(saved["database"])) as db:
            self.assertEqual(db.execute("SELECT count(*) FROM runs").fetchone()[0], 1)
            self.assertEqual(db.execute("SELECT count(*) FROM reference_bgcs").fetchone()[0], 2)
            self.assertEqual(db.execute("SELECT count(*) FROM bgc_compounds").fetchone()[0], 4)
            links = db.execute('''SELECT r.region_id, r.accession, r.relation, c.inchikey
                FROM region_bgcs r JOIN bgc_compounds b ON r.run_id=b.run_id AND r.accession=b.accession
                JOIN compounds c ON c.compound_id=b.compound_id WHERE c.inchikey IS NOT NULL''').fetchall()
            self.assertEqual(links, [("r1", "BGC0000001", "same_gcf", KEY), ("r1", "BGC0000002", "same_gcf", KEY)])
            self.assertEqual(json.loads(db.execute("SELECT provenance FROM runs").fetchone()[0]), first["provenance"])
        self.assertEqual(json.loads(Path(saved["evidence"]).read_text())["compounds"], first["compounds"])

    def test_multiple_identifiers_require_consistent_keys(self):
        self.records["BGC0000001"]["compounds"][0]["databaseIds"] = ["pubchem:2244", "pubchem:2519"]
        write_metadata(self.archive, self.records)
        with patch("compound_context.get", side_effect=[self.response(), self.response(2519, OTHER_KEY)]):
            result = self.enrich()
        self.assertIsNone(result["compounds"][result["bgc_compounds"][0][1]])
        (self.root / "pubchem/2519.json").unlink()
        with patch("compound_context.get", return_value=self.response(2519, KEY)):
            consistent = self.enrich()
        self.assertEqual(consistent["bgc_compounds"][0][1], KEY)

    def test_service_failure_is_recorded_without_discarding_bgcs_and_can_be_retried(self):
        with patch("compound_context.get", side_effect=requests.Timeout("service unavailable")):
            result = self.enrich()
        self.assertEqual(len(result["reference_bgcs"]), 2)
        self.assertTrue(all(key is None for key in result["compounds"].values()))
        self.assertIn("service unavailable", result["provenance"]["compounds"]["pubchem"]["2244"]["error"])
        self.assertFalse((self.root / "pubchem/2244.json").exists())
        save_results(result, self.root)
        with patch("compound_context.get", return_value=self.response()):
            recovered = self.enrich()
        saved = save_results(recovered, self.root)
        with closing(sqlite3.connect(saved["database"])) as db:
            self.assertEqual(db.execute("SELECT count(*) FROM runs").fetchone()[0], 2)

    def test_bad_cid_or_key_is_not_persisted_as_a_resolved_compound(self):
        for response in [self.response(cid=999), self.response(key="invalid")]:
            with self.subTest(response=response), patch("compound_context.get", return_value=response):
                result = self.enrich()
            self.assertTrue(all(key is None for key in result["compounds"].values()))
            self.assertIn("error", result["provenance"]["compounds"]["pubchem"]["2244"])

    def test_no_assignment_makes_no_compound_claims_or_requests(self):
        self.evidence["regions"][0].update(gcf_id=None, gcf_status="no_features")
        with patch("compound_context.get") as get:
            result = self.enrich()
        get.assert_not_called()
        self.assertEqual(result["reference_bgcs"], {})
        self.assertEqual(result["compounds"], {})
        self.assertEqual(result["region_bgcs"], [])

    def test_metadata_release_and_accession_must_match_reference(self):
        write_metadata(self.archive, self.records, release="3.1")
        with self.assertRaises(KeyError):
            self.enrich()
        self.records["BGC0000001"]["accession"] = "BGC0000002"
        write_metadata(self.archive, self.records)
        with self.assertRaisesRegex(ValueError, "accession mismatch"):
            self.enrich()

    def test_stage_and_ledger_stay_within_seventy_five_lines(self):
        root = Path(__file__).resolve().parents[1]
        for name in ("compound_context.py", "save_results.py"):
            self.assertLessEqual(len((root / name).read_text().splitlines()), 75)
