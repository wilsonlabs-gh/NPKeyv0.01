"""Stage 3: query BGCs against a frozen MIBiG-derived BiG-SLiCE reference."""

import hashlib
import json
import math
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
from .io_utils import checksum, execute, rows, write_json


def assign_gcfs(evidence, reference_dir, output_dir, *, cutoff=None, ranks=3, executable="bigslice"):
    reference = Path(reference_dir).expanduser().resolve()
    model = json.loads((reference / "reference.json").read_text(encoding="utf-8"))
    cutoff = model["threshold"] if cutoff is None else cutoff
    if not math.isfinite(cutoff) or cutoff < 0 or not isinstance(ranks, int) or ranks < 1:
        raise ValueError("Use a finite non-negative cutoff and a positive integer ranks.")
    source_db = reference / "engine/result/data.db"
    if checksum(source_db) != model["model_id"]:
        raise ValueError("Reference database checksum changed.")
    source_cache = reference / "engine/result/cache" / model["cache_file"]
    if checksum(source_cache) != model["cache_checksum"]:
        raise ValueError("Reference model cache checksum changed.")
    ranks = min(ranks, len(model["families"]))
    if not ranks:
        raise ValueError("Reference has no families.")
    regions = evidence["regions"]
    for region in regions:
        if checksum(region["bgc_path"]) != region["bgc_checksum"]:
            raise ValueError(f"BGC file changed: {region['bgc_path']}")
    inputs = [(region["region_id"], region["bgc_checksum"]) for region in regions]
    identity = [model["model_id"], cutoff, ranks, sorted(inputs)]
    query_id = hashlib.sha256(json.dumps(identity).encode()).hexdigest()
    work = Path(output_dir).expanduser().resolve() / "gcf" / query_id
    work.mkdir(parents=True, exist_ok=True)
    query_db = work / "query.db"
    if regions and not query_db.exists():
        execute([executable, "--version"], work / "version.txt")
        if (work / "version.txt").read_text().strip() != model["engine_version"]:
            raise ValueError("BiG-SLiCE version/HMM setup differs from the reference build.")
        with TemporaryDirectory(prefix="query-", dir=work) as temporary:
            stage = Path(temporary)
            (stage / "engine/result/cache").mkdir(parents=True)
            (stage / "input").mkdir()
            shutil.copyfile(source_db, stage / "engine/result/data.db")
            shutil.copyfile(source_cache, stage / "engine/result/cache" / model["cache_file"])
            for region in regions:
                target = stage / "input" / f"{region['region_id']}.region001.gbk"
                shutil.copyfile(region["bgc_path"], target)
            execute([executable, "--query", stage / "input", "--query_name", "query",
                     "--run_id", str(model["run_id"]), "--n_ranks", str(ranks), "-t", "1",
                     "--program_db_folder", model["hmm_dir"], stage / "engine"], work / "query.log")
            report = rows(stage / "engine/reports/reports.db",
                          "SELECT id FROM reports WHERE type='query' AND name='query'")
            if len(report) != 1:
                raise ValueError("Expected one completed BiG-SLiCE query report.")
            raw = stage / f"engine/reports/{report[0]['id']}/data.db"
            _read_matches(raw, regions, model, cutoff)
            shutil.copyfile(raw, work / "query.tmp")
            (work / "query.tmp").replace(query_db)
    evidence["gcf_matches"] = _read_matches(query_db, regions, model, cutoff) if regions else []
    evidence["provenance"]["gcf"] = dict(model, cutoff=cutoff, ranks=ranks,
                                         query_db=str(query_db) if regions else None)
    write_json(work / "evidence.json", evidence)
    return evidence


def _read_matches(query_db, regions, model, cutoff):
    queries = rows(query_db, "SELECT b.id, b.orig_filename, EXISTS(SELECT 1 FROM bgc_features f "
                   "WHERE f.bgc_id=b.id AND f.value>0) AS has_features FROM bgc b")
    mapping = {q["id"]: q["orig_filename"].split(".region")[0] for q in queries}
    by_id = {region["region_id"]: region for region in regions}
    if len(queries) != len(regions) or set(mapping.values()) != set(by_id):
        raise ValueError("BiG-SLiCE did not return every input BGC exactly once.")
    featured = {q["id"] for q in queries if q["has_features"]}
    for region in regions:
        region.update(gcf_id=None, gcf_status="unassigned")
    matches = []
    first_hits = set()
    for hit in rows(query_db, "SELECT * FROM gcf_membership ORDER BY bgc_id, rank"):
        region_id = mapping[hit["bgc_id"]]
        family, distance = str(hit["gcf_id"]), hit["membership_value"]
        if family not in model["families"] or not math.isfinite(distance) or distance < 0:
            raise ValueError("Invalid GCF identifier or distance in BiG-SLiCE output.")
        accepted = hit["bgc_id"] in featured and distance <= cutoff
        if hit["rank"] == 0:
            first_hits.add(region_id)
            status = "assigned" if accepted else "unassigned"
            if hit["bgc_id"] not in featured:
                status = "no_features"
            by_id[region_id].update(gcf_id=hit["gcf_id"] if accepted else None, gcf_status=status)
        matches.append({"region_id": region_id, "gcf_id": hit["gcf_id"], "rank": hit["rank"] + 1,
                        "distance": distance, "within_cutoff": accepted})
    if first_hits != set(by_id):
        raise ValueError("Missing nearest-family results in BiG-SLiCE output.")
    return sorted(matches, key=lambda hit: (hit["region_id"], hit["rank"]))
