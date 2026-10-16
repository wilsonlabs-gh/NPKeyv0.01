"""Stage 5: persist BGC/compound nodes, relationships and provenance."""

from contextlib import closing
import csv
import hashlib
import json
from pathlib import Path
import sqlite3
from io_utils import write_json


def save_results(evidence, output_dir):
    run_id = hashlib.sha256(json.dumps(evidence, sort_keys=True).encode()).hexdigest()
    root = Path(output_dir).expanduser().resolve()
    folder = root / "results" / run_id
    folder.mkdir(parents=True, exist_ok=True)
    model = evidence["provenance"]["gcf"]
    database = root / "ledger.sqlite"
    with closing(sqlite3.connect(database)) as db, db:
        db.executescript("""
            CREATE TABLE IF NOT EXISTS runs (
                run_id TEXT PRIMARY KEY, provenance TEXT, issues TEXT);
            CREATE TABLE IF NOT EXISTS regions (
                run_id TEXT, region_id TEXT, sequence_id TEXT, start INTEGER, end INTEGER,
                gcf_id INTEGER, status TEXT, details TEXT, PRIMARY KEY(run_id, region_id));
            CREATE TABLE IF NOT EXISTS families (
                model_id TEXT, gcf_id INTEGER, mibig_accessions TEXT,
                PRIMARY KEY(model_id, gcf_id));
            CREATE TABLE IF NOT EXISTS gcf_matches (
                run_id TEXT, region_id TEXT, model_id TEXT, gcf_id INTEGER, rank INTEGER,
                distance REAL, within_cutoff INTEGER, PRIMARY KEY(run_id, region_id, rank));
            CREATE TABLE IF NOT EXISTS reference_bgcs (
                run_id TEXT, accession TEXT, metadata TEXT, PRIMARY KEY(run_id, accession));
            CREATE TABLE IF NOT EXISTS compounds (
                compound_id TEXT PRIMARY KEY, inchikey TEXT);
            CREATE TABLE IF NOT EXISTS bgc_compounds (
                run_id TEXT, accession TEXT, compound_id TEXT, source_index INTEGER, evidence TEXT,
                PRIMARY KEY(run_id, accession, source_index));
            CREATE TABLE IF NOT EXISTS region_bgcs (
                run_id TEXT, region_id TEXT, accession TEXT, relation TEXT, PRIMARY KEY(run_id, region_id, accession));
        """)
        db.execute("INSERT OR REPLACE INTO runs VALUES (?,?,?)",
                   (run_id, json.dumps(evidence["provenance"]), json.dumps(evidence["issues"])))
        for region in evidence["regions"]:
            db.execute("INSERT OR REPLACE INTO regions VALUES (?,?,?,?,?,?,?,?)",
                       (run_id, region["region_id"], region["sequence_id"], region["start"],
                        region["end"], region["gcf_id"], region["gcf_status"], json.dumps(region)))
        db.executemany("INSERT OR REPLACE INTO families VALUES (?,?,?)",
                       [(model["model_id"], int(family), json.dumps(members)) for family, members in model["families"].items()])
        db.executemany("INSERT OR REPLACE INTO gcf_matches VALUES (?,?,?,?,?,?,?)",
                       [(run_id, hit["region_id"], model["model_id"], hit["gcf_id"], hit["rank"],
                         hit["distance"], hit["within_cutoff"]) for hit in evidence["gcf_matches"]])
        db.executemany("INSERT OR REPLACE INTO reference_bgcs VALUES (?,?,?)",
                       [(run_id, accession, json.dumps(record)) for accession, record in evidence.get("reference_bgcs", {}).items()])
        db.executemany("INSERT OR REPLACE INTO compounds VALUES (?,?)", evidence.get("compounds", {}).items())
        db.executemany("INSERT OR REPLACE INTO bgc_compounds VALUES (?,?,?,?,?)",
                       [(run_id, accession, compound, index, json.dumps(detail)) for accession, compound, index, detail in evidence.get("bgc_compounds", [])])
        db.executemany("INSERT OR REPLACE INTO region_bgcs VALUES (?,?,?,?)",
                       [(run_id, region, accession, "same_gcf") for region, accession in evidence.get("region_bgcs", [])])
    write_json(folder / "evidence.json", evidence)
    nearest = {hit["region_id"]: hit for hit in evidence["gcf_matches"] if hit["rank"] == 1}
    fields = ["region_id", "sequence_id", "start", "end", "gcf_id", "gcf_status",
              "nearest_gcf_id", "distance", "cutoff", "model_id", "mibig_accessions"]
    with (folder / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for region in evidence["regions"]:
            row = {key: region[key] for key in fields[:6]}
            hit = nearest.get(region["region_id"], {})
            row.update(nearest_gcf_id=hit.get("gcf_id"), distance=hit.get("distance"),
                       cutoff=model["cutoff"], model_id=model["model_id"],
                       mibig_accessions=";".join(model["families"].get(str(region["gcf_id"]), [])))
            writer.writerow(row)
    return {"run_id": run_id, "status": "complete", "database": str(database),
            "evidence": str(folder / "evidence.json"), "summary": str(folder / "summary.csv")}
