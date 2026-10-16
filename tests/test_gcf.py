"""Synthetic SQLite contracts from BiG-SLiCE e622dabf; external HMM work is mocked."""

from contextlib import closing
from copy import deepcopy
from io import BytesIO
import json
from pathlib import Path
import sqlite3
import subprocess
import tarfile
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from antismash import run_antismash
from gcf_assignment import assign_gcfs
from gcf_reference import prepare_reference
from io_utils import checksum
from main import run
from save_results import save_results
from test_antismash import result_data, region_files
from test_ingest import genbank, record
from test_compounds import write_metadata

VERSION = "BiG-SLiCE version 2.0.0\nHMM databases version: synthetic"


def model_database(path, filenames, threshold=0.4, status=7):
    path.parent.mkdir(parents=True, exist_ok=True)
    with closing(sqlite3.connect(path)) as db, db:
        db.executescript('''
            CREATE TABLE run (id INTEGER, status INTEGER, hmm_db_id INTEGER);
            CREATE TABLE hmm_db (id INTEGER, md5_biosyn_pfam TEXT, md5_sub_pfam TEXT);
            CREATE TABLE clustering (id INTEGER, run_id INTEGER, threshold REAL, random_seed INTEGER);
            CREATE TABLE bgc (id INTEGER, orig_filename TEXT);
            CREATE TABLE gcf (id INTEGER, clustering_id INTEGER);
            CREATE TABLE gcf_membership (gcf_id INTEGER, bgc_id INTEGER, rank INTEGER, membership_value REAL);
        ''')
        db.execute("INSERT INTO run VALUES (1,?,7)", (status,))
        db.execute("INSERT INTO hmm_db VALUES (7,'synthetic-pfam','synthetic-subpfam')")
        db.execute("INSERT INTO clustering VALUES (1,1,?,123)", (threshold,))
        for number, filename in enumerate(filenames, 1):
            db.execute("INSERT INTO bgc VALUES (?,?)", (number, filename))
            db.execute("INSERT INTO gcf VALUES (?,1)", (number,))
            db.execute("INSERT INTO gcf_membership VALUES (?,?,0,0.0)", (number, number))


class GcfTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.reference = self.root / "reference"
        self.archive = self.root / "release.tar.bz2"
        self.hmm = self.root / "hmms"
        self.hmm.mkdir()
        self.output = self.root / "output"
        self.distances = [0.2, 0.6]
        self.has_features = True
        self.status = 7
        self.release()

    def release(self, prefix="mibig_antismash_4.0_gbk_as8b1"):
        with tarfile.open(self.archive, "w:bz2") as archive:
            for accession in ("BGC0000001", "BGC0000002"):
                content = genbank(record(accession))
                info = tarfile.TarInfo(f"{prefix}/{accession}.gbk")
                info.size = len(content)
                archive.addfile(info, BytesIO(content))

    def build_command(self, command, log):
        if "--version" in command:
            Path(log).write_text(VERSION)
            return
        self.assertEqual(command[command.index("--n_ranks") + 1], "1")
        self.assertNotIn("--complete", command)  # MIBiG subregions can be labelled edge=True.
        self.assertIn("--early_dumping", command)
        inputs = self.reference / "input/mibig/records"
        files = sorted(path.name for path in inputs.glob("*.gbk"))
        database = Path(command[-1]) / "result/data.db"
        if "--resume" in command:
            with closing(sqlite3.connect(database)) as db, db:
                db.execute("UPDATE run SET status=?", (self.status,))
        else:
            model_database(database, files, status=self.status)
        cache = Path(command[-1]) / "result/cache"
        cache.mkdir(exist_ok=True)
        (cache / "gcf_models_1.pkl").write_bytes(b"synthetic centroid cache")
        Path(log).write_text("Synthetic reference build")

    def build(self):
        with patch("gcf_reference.execute", side_effect=self.build_command):
            return prepare_reference(self.archive, self.reference, self.hmm)

    def evidence(self):
        source = self.root / "genome.gbk"
        source.write_bytes(genbank(record("TEST.1", "ATG" * 40)))
        saved = self.root / "antismash"
        saved.mkdir(exist_ok=True)
        data = result_data()
        (saved / "genome.json").write_text(json.dumps(data))
        for name, contents in region_files(data).items():
            (saved / name).write_text(contents)
        from ingest import ingest_genome
        genome = ingest_genome(source, self.output)
        return run_antismash(genome, self.output, saved_results=saved)

    def query_command(self, command, log):
        if "--version" in command:
            Path(log).write_text(VERSION)
            return
        self.assertIn("--run_id", command)
        self.assertEqual(command[command.index("--run_id") + 1], "1")
        inputs = Path(command[command.index("--query") + 1])
        engine = Path(command[-1])
        self.assertNotEqual(engine, self.reference / "engine")
        self.assertEqual((engine / "result/cache/gcf_models_1.pkl").read_bytes(), b"synthetic centroid cache")
        reports = engine / "reports"
        (reports / "1").mkdir(parents=True)
        with closing(sqlite3.connect(reports / "reports.db")) as db, db:
            db.execute("CREATE TABLE reports (id INTEGER, type TEXT, name TEXT)")
            db.execute("INSERT INTO reports VALUES (1,'query','query')")
        with closing(sqlite3.connect(reports / "1/data.db")) as db, db:
            db.executescript('''
                CREATE TABLE bgc (id INTEGER, orig_filename TEXT);
                CREATE TABLE bgc_features (bgc_id INTEGER, value INTEGER);
                CREATE TABLE gcf_membership (bgc_id INTEGER, gcf_id INTEGER,
                                            rank INTEGER, membership_value REAL);
            ''')
            for number, path in enumerate(sorted(inputs.glob("*.gbk")), 1):
                db.execute("INSERT INTO bgc VALUES (?,?)", (number, path.name))
                if self.has_features:
                    db.execute("INSERT INTO bgc_features VALUES (?,255)", (number,))
                for rank, distance in enumerate(self.distances):
                    db.execute("INSERT INTO gcf_membership VALUES (?,?,?,?)",
                               (number, rank + 1, rank, distance))
        Path(log).write_text("Synthetic query")

    def assign(self, evidence, **options):
        with patch("gcf_assignment.execute", side_effect=self.query_command):
            return assign_gcfs(evidence, self.reference, self.output, **options)

    def test_reference_records_families_version_hmm_and_exact_source(self):
        self.assertEqual(self.build(), self.reference)
        manifest = json.loads((self.reference / "reference.json").read_text())
        self.assertEqual(manifest["mibig_release"], "4.0")
        self.assertEqual(manifest["archive_checksum"], checksum(self.archive))
        self.assertEqual(manifest["run_id"], 1)
        self.assertEqual(manifest["clustering"]["id"], 7)
        self.assertEqual(manifest["families"], {"1": ["BGC0000001"], "2": ["BGC0000002"]})
        for path in (self.reference / "input/mibig/records").glob("*.gbk"):
            self.assertEqual(path.read_bytes(), genbank(record(path.stem)))

    def test_reference_reuse_and_changed_inputs(self):
        self.build()
        with patch("gcf_reference.execute", side_effect=self.build_command) as engine:
            self.assertEqual(prepare_reference(self.archive, self.reference, self.hmm), self.reference)
        self.assertEqual(engine.call_count, 1)  # Version check only.
        with patch("gcf_reference.execute", side_effect=self.build_command):
            with self.assertRaisesRegex(ValueError, "inputs changed"):
                prepare_reference(self.archive, self.reference, self.hmm, threshold=0.5)

    def test_incomplete_reference_is_not_published(self):
        self.status = 5
        with self.assertRaisesRegex(ValueError, "completed"):
            self.build()
        self.assertFalse((self.reference / "reference.json").exists())
        self.status = 7
        with patch("gcf_reference.execute", side_effect=self.build_command) as engine:
            prepare_reference(self.archive, self.reference, self.hmm)
        self.assertIn("--resume", engine.call_args.args[0])
        self.assertTrue((self.reference / "reference.json").exists())

    def test_centroid_cache_integrity_is_checked(self):
        self.build()
        cache = self.reference / "engine/result/cache/gcf_models_1.pkl"
        cache.write_bytes(b"changed centroids")
        with self.assertRaisesRegex(ValueError, "cache checksum"):
            self.assign(self.evidence())
        with self.assertRaisesRegex(ValueError, "cache checksum"):
            self.build()

    def test_multiple_bgcs_are_linked_by_id_without_duplicate_family_records(self):
        self.build()
        evidence = self.evidence()
        second = deepcopy(evidence["regions"][0])
        second["region_id"] = "independent-query"
        second["sequence_id"] = "OTHER.2"
        evidence["regions"].append(second)
        result = self.assign(evidence)
        self.assertEqual({hit["region_id"] for hit in result["gcf_matches"]},
                         {region["region_id"] for region in evidence["regions"]})
        self.assertEqual(len(result["gcf_matches"]), 4)
        saved = save_results(result, self.output)
        with closing(sqlite3.connect(saved["database"])) as db:
            self.assertEqual(db.execute("SELECT count(*) FROM regions").fetchone()[0], 2)
            self.assertEqual(db.execute("SELECT count(*) FROM families").fetchone()[0], 2)

    def test_raw_mibig_export_is_rejected(self):
        self.release(prefix="mibig_gbk_4.0")
        with self.assertRaisesRegex(ValueError, "antiSMASH-processed"):
            self.build()

    def test_handoff_assigns_by_distance_keeps_candidates_and_reference_frozen(self):
        self.build()
        evidence = self.evidence()
        original = deepcopy(evidence)
        database = self.reference / "engine/result/data.db"
        before = checksum(database)
        result = self.assign(evidence)
        self.assertEqual(checksum(database), before)
        region = result["regions"][0]
        self.assertEqual(region["gcf_status"], "assigned")
        self.assertEqual(region["gcf_id"], 1)
        for key, value in original["regions"][0].items():
            self.assertEqual(region[key], value)
        self.assertEqual([hit["rank"] for hit in result["gcf_matches"]], [1, 2])
        self.assertEqual([hit["within_cutoff"] for hit in result["gcf_matches"]], [True, False])
        self.assertEqual(result["provenance"]["gcf"]["ranks"], 2)
        first = deepcopy(result)
        with patch("gcf_assignment.execute") as engine:
            self.assertEqual(assign_gcfs(result, self.reference, self.output), first)
        engine.assert_not_called()

    def test_unassigned_and_no_features_keep_all_bgc_candidates(self):
        self.build()
        self.distances = [0.41, 0.9]
        result = self.assign(self.evidence())
        self.assertEqual(result["regions"][0]["gcf_status"], "unassigned")
        self.assertIsNone(result["regions"][0]["gcf_id"])
        self.assertEqual(len(result["gcf_matches"]), 2)
        self.has_features = False
        self.distances = [0.0, 0.0]
        result = self.assign(self.evidence(), cutoff=0.5)
        self.assertEqual(result["regions"][0]["gcf_status"], "no_features")
        self.assertTrue(all(not hit["within_cutoff"] for hit in result["gcf_matches"]))

    def test_exact_cutoff_is_accepted_and_changes_have_separate_cache(self):
        self.build()
        self.distances = [0.4, 0.8]
        first = deepcopy(self.assign(self.evidence()))
        second = self.assign(self.evidence(), cutoff=0.39)
        self.assertEqual(first["regions"][0]["gcf_status"], "assigned")
        self.assertEqual(second["regions"][0]["gcf_status"], "unassigned")
        self.assertNotEqual(first["provenance"]["gcf"]["query_db"], second["provenance"]["gcf"]["query_db"])

    def test_no_bgcs_is_valid_without_invoking_engine(self):
        self.build()
        evidence = self.evidence()
        evidence["regions"] = []
        with patch("gcf_assignment.execute") as engine:
            result = assign_gcfs(evidence, self.reference, self.output)
        engine.assert_not_called()
        self.assertEqual(result["gcf_matches"], [])
        self.assertIsNone(result["provenance"]["gcf"]["query_db"])
        self.assertEqual(save_results(result, self.output)["status"], "complete")

    def test_changed_bgc_or_model_is_rejected(self):
        self.build()
        evidence = self.evidence()
        with Path(evidence["regions"][0]["bgc_path"]).open("a") as handle:
            handle.write("\n")
        with self.assertRaisesRegex(ValueError, "BGC file changed"):
            self.assign(evidence)
        evidence = self.evidence()
        with (self.reference / "engine/result/data.db").open("ab") as handle:
            handle.write(b"modified")
        with self.assertRaisesRegex(ValueError, "database checksum"):
            self.assign(evidence)

    def test_failed_command_keeps_prediction_and_has_no_completed_query(self):
        self.build()
        evidence = self.evidence()
        original = deepcopy(evidence)
        def fail(command, log):
            if "--version" in command:
                Path(log).write_text(VERSION)
            else:
                raise subprocess.CalledProcessError(1, command)
        with patch("gcf_assignment.execute", side_effect=fail):
            with self.assertRaises(subprocess.CalledProcessError):
                assign_gcfs(evidence, self.reference, self.output)
        self.assertEqual(evidence, original)
        self.assertEqual(list(self.output.rglob("query.db")), [])
        self.assertTrue(list(self.output.glob("antismash/*/evidence.json")))
        self.assertEqual(self.assign(evidence)["regions"][0]["gcf_status"], "assigned")

    def test_bad_memberships_and_missing_query_rows_are_explicit(self):
        self.build()
        result = self.assign(self.evidence())
        path = Path(result["provenance"]["gcf"]["query_db"])
        with closing(sqlite3.connect(path)) as db, db:
            db.execute("UPDATE gcf_membership SET gcf_id=999 WHERE rank=0")
        with self.assertRaisesRegex(ValueError, "Invalid GCF"):
            self.assign(self.evidence())
        with closing(sqlite3.connect(path)) as db, db:
            db.execute("DELETE FROM bgc")
        with self.assertRaisesRegex(ValueError, "every input BGC"):
            self.assign(self.evidence())

    def test_full_controller_saves_reference_compounds_repeatably(self):
        self.build()
        self.evidence()
        metadata = self.root / "metadata.tar.gz"
        write_metadata(metadata, {"BGC0000001": {"accession": "BGC0000001", "version": 1,
                       "compounds": [{"name": "unresolved reference compound", "databaseIds": []}]}})
        with patch("gcf_assignment.execute", side_effect=self.query_command):
            first = run(self.root / "genome.gbk", self.output, self.reference,
                        saved_results=self.root / "antismash", mibig_metadata=metadata)
            second = run(self.root / "genome.gbk", self.output, self.reference,
                         saved_results=self.root / "antismash", mibig_metadata=metadata)
        self.assertEqual(first, second)
        with closing(sqlite3.connect(first["database"])) as db:
            self.assertEqual(db.execute("SELECT count(*) FROM runs").fetchone()[0], 1)
            self.assertEqual(db.execute("SELECT count(*) FROM regions").fetchone()[0], 1)
            self.assertEqual(db.execute("SELECT count(*) FROM gcf_matches").fetchone()[0], 2)
            self.assertEqual(db.execute("SELECT count(*) FROM bgc_compounds").fetchone()[0], 1)
        saved = json.loads(Path(first["evidence"]).read_text())
        self.assertEqual(len(saved["compounds"]), 1)
        self.assertIn("BGC0000001", saved["reference_bgcs"])
        self.assertIn("BGC0000001", Path(first["summary"]).read_text())

    def test_invalid_settings_are_not_silently_coerced(self):
        self.build()
        for options in ({"cutoff": float("nan")}, {"cutoff": -1}, {"ranks": 0}, {"ranks": 1.5}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.assign(self.evidence(), **options)

    def test_production_files_stay_under_one_hundred_physical_lines(self):
        root = Path(__file__).resolve().parents[1]
        for path in root.glob("*.py"):
            with self.subTest(file=path.name):
                self.assertLess(len(path.read_text(encoding="utf-8").splitlines()), 100)
