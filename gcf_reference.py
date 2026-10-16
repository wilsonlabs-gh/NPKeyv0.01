"""Build BiG-SLiCE families from an official antiSMASH-processed MIBiG archive."""

import json
import math
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
from io_utils import checksum, execute, rows, write_json


def prepare_reference(mibig_archive, output_dir, hmm_dir, *, threshold=0.4, executable="bigslice"):
    if not math.isfinite(threshold) or threshold <= 0:
        raise ValueError("Clustering threshold must be finite and positive.")
    root = Path(output_dir).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    archive_path = Path(mibig_archive).expanduser().resolve()
    hmm_dir = Path(hmm_dir).expanduser().resolve()
    execute([executable, "--version"], root / "version.txt")
    request = {"archive_checksum": checksum(archive_path), "threshold": threshold,
               "hmm_dir": str(hmm_dir),
               "engine_version": (root / "version.txt").read_text(encoding="utf-8").strip()}
    request_path = root / "request.json"
    if request_path.exists() and json.loads(request_path.read_text()) != request:
        raise ValueError("Reference inputs changed; use a new reference directory.")
    write_json(request_path, request)
    manifest_path = root / "reference.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if checksum(root / "engine/result/data.db") != manifest["model_id"]:
            raise ValueError("Reference database checksum changed.")
        if checksum(root / "engine/result/cache" / manifest["cache_file"]) != manifest["cache_checksum"]:
            raise ValueError("Reference model cache checksum changed.")
        return root
    inputs = root / "input"
    records = inputs / "mibig" / "records"
    records.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path) as archive:
        members = [m for m in archive.getmembers() if m.isfile()
                   and re.fullmatch(r"BGC[0-9]+\.gbk", PurePosixPath(m.name).name)]
        folders = {PurePosixPath(m.name).parent.name for m in members}
        release_match = re.fullmatch(r"mibig_antismash_([0-9.]+)_gbk.*", next(iter(folders), ""))
        if len(folders) != 1 or not release_match:
            raise ValueError("Expected an official antiSMASH-processed MIBiG GenBank archive.")
        release = release_match[1]
        names = [PurePosixPath(m.name).name for m in members]
        if len(set(names)) != len(names):
            raise ValueError("Duplicate MIBiG accession in archive.")
        for member, name in zip(members, names):
            with archive.extractfile(member) as source, (records / name).open("wb") as target:
                shutil.copyfileobj(source, target)
    (inputs / "datasets.tsv").write_text(f"mibig\tmibig\t\tMIBiG {release}\n", encoding="utf-8")
    engine = root / "engine"
    options = ["--resume", "--run_id", "1"] if engine.exists() else ["-i", inputs]
    execute([executable, *options, "--threshold", str(threshold), "--n_ranks", "1",
             "--early_dumping", "-t", "1", "--program_db_folder", hmm_dir, engine], root / "build.log")
    db = engine / "result/data.db"
    runs = rows(db, "SELECT r.id AS run_id, r.status, c.id AS clustering_id, c.threshold, c.random_seed, h.* FROM run r "
                "JOIN clustering c ON c.run_id=r.id JOIN hmm_db h ON h.id=r.hmm_db_id")
    if len(runs) != 1 or runs[0]["status"] != 7:
        raise ValueError("Expected one completed BiG-SLiCE reference run.")
    parsed = rows(db, "SELECT orig_filename FROM bgc")
    if sorted(row["orig_filename"] for row in parsed) != sorted(names):
        raise ValueError("BiG-SLiCE did not parse every MIBiG record; inspect build.log.")
    families = {row["id"]: [] for row in rows(db, "SELECT id FROM gcf")}
    memberships = rows(db, "SELECT m.gcf_id, b.orig_filename FROM gcf_membership m "
                       "JOIN bgc b ON b.id=m.bgc_id WHERE m.rank=0 AND m.membership_value<=?",
                       (threshold,))
    for member in memberships:
        families[member["gcf_id"]].append(Path(member["orig_filename"]).stem)
    if not families:
        raise ValueError("BiG-SLiCE produced no reference families.")
    cache_file = f"gcf_models_{runs[0]['clustering_id']}.pkl"
    manifest = dict(request, model_id=checksum(db), mibig_release=release, cache_file=cache_file,
                    cache_checksum=checksum(engine / "result/cache" / cache_file),
                    source_archive=str(archive_path), run_id=runs[0]["run_id"],
                    clustering=runs[0],
                    families={str(key): sorted(value) for key, value in families.items()})
    write_json(manifest_path, manifest)
    return root
