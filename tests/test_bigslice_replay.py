"""Replay published BiG-SLiCE 2 results through our reader and saved outputs."""

from contextlib import closing
from copy import deepcopy
import csv
import gzip
import json
import math
from pathlib import Path
import sqlite3
from tempfile import TemporaryDirectory
import unittest

from npkey.gcf_assignment import _read_matches
from npkey.io_utils import checksum
from npkey.save_results import save_results

FIXTURE = Path(__file__).parent / "fixtures/bigslice_v2"


def example(root):
    metadata = json.loads((FIXTURE / "source.json").read_text())
    database = root / "query.db"
    database.write_bytes(gzip.decompress((FIXTURE / "query.db.gz").read_bytes()))
    model = metadata["model"]
    regions = [{"region_id": row["orig_filename"].split(".region")[0],
                "sequence_id": row["orig_filename"], "start": 0, "end": row["length_nt"]}
               for row in metadata["bgcs"]]
    evidence = {"regions": regions, "gcf_matches": [], "issues": [],
                "provenance": {"gcf": model, "fixture_source": metadata["source"],
                               "purpose": "Published BGC output replay; coordinates are BGC-local."}}
    return database, evidence, metadata


class PublishedOutputTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.database, self.evidence, self.metadata = example(self.root)
        self.model = self.evidence["provenance"]["gcf"]

    def read(self, cutoff=None):
        return _read_matches(self.database, self.evidence["regions"], self.model,
                             self.model["cutoff"] if cutoff is None else cutoff)

    def test_exact_query_schema_preserves_published_values(self):
        self.assertEqual(checksum(self.database), self.metadata["query_sha256"])
        with closing(sqlite3.connect(self.database)) as db:
            self.assertEqual(db.execute("SELECT * FROM schema").fetchall(), [("1.0.1", "1.0.1")])
            columns = {row[1]: row[2] for row in db.execute("PRAGMA table_info(gcf_membership)")}
            self.assertEqual(columns["membership_value"], "INTEGER")
            actual = db.execute("SELECT bgc_id,gcf_id,rank,membership_value FROM gcf_membership "
                                "ORDER BY bgc_id,rank").fetchall()
        expected = sorted((bgc["id"], hit["gcf_id"], hit["rank"], hit["membership_value"])
                          for bgc in self.metadata["bgcs"] for hit in bgc["memberships"])
        self.assertEqual(actual, expected)
        self.assertTrue(any(0 < row[3] < 1 for row in actual))

    def test_published_ids_ranks_and_distances_survive_all_exports_and_replay(self):
        before = checksum(self.database)
        self.evidence["regions"].reverse()
        self.evidence["gcf_matches"] = self.read()
        expected = []
        for bgc in self.metadata["bgcs"]:
            for hit in bgc["memberships"]:
                expected.append({"region_id": bgc["orig_filename"].split(".region")[0],
                                 "gcf_id": hit["gcf_id"], "rank": hit["rank"] + 1,
                                 "distance": hit["membership_value"],
                                 "within_cutoff": bgc["has_features"] and
                                 hit["membership_value"] <= self.model["cutoff"]})
        self.assertEqual(self.evidence["gcf_matches"], sorted(expected, key=lambda x: (x["region_id"], x["rank"])))
        self.assertEqual(sorted(row["gcf_status"] for row in self.evidence["regions"]),
                         ["assigned", "assigned", "no_features", "unassigned"])
        saved = save_results(self.evidence, self.root)
        self.assertEqual(save_results(deepcopy(self.evidence), self.root), saved)
        self.assertEqual(json.loads(Path(saved["evidence"]).read_text()), self.evidence)
        with closing(sqlite3.connect(saved["database"])) as db:
            actual = db.execute("SELECT region_id,gcf_id,rank,distance,within_cutoff "
                                "FROM gcf_matches ORDER BY region_id,rank").fetchall()
            self.assertEqual(actual, [tuple(hit.values()) for hit in self.evidence["gcf_matches"]])
            self.assertEqual(db.execute("SELECT count(*) FROM runs").fetchone()[0], 1)
        with Path(saved["summary"]).open(newline="") as handle:
            summary = {row["region_id"]: row for row in csv.DictReader(handle)}
        for hit in self.evidence["gcf_matches"]:
            if hit["rank"] == 1:
                self.assertEqual(float(summary[hit["region_id"]]["distance"]), hit["distance"])
                self.assertEqual(int(summary[hit["region_id"]]["nearest_gcf_id"]), hit["gcf_id"])
        self.assertEqual(checksum(self.database), before)

    def test_cutoff_boundary_on_a_real_fractional_distance(self):
        bgc = next(row for row in self.metadata["bgcs"] if row["has_features"] and
                   0 < row["memberships"][0]["membership_value"] < 1)
        distance = bgc["memberships"][0]["membership_value"]
        identity = bgc["orig_filename"].split(".region")[0]
        self.read(cutoff=distance)
        region = next(row for row in self.evidence["regions"] if row["region_id"] == identity)
        self.assertEqual(region["gcf_status"], "assigned")
        self.read(cutoff=math.nextafter(distance, 0))
        self.assertEqual(region["gcf_status"], "unassigned")
