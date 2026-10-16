"""Parse antiSMASH 8 regions and verify the BGC files passed to BiG-SLiCE."""

import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from io_utils import checksum


def parse_results(results_dir, genome):
    results_dir = Path(results_dir).expanduser().resolve()
    candidates = []
    for path in results_dir.rglob("*.json"):
        data = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(data, dict) and "records" in data and "version" in data:
            candidates.append((path, data))
    if len(candidates) != 1:
        raise ValueError("Expected one antiSMASH results JSON file.")
    result_path, data = candidates[0]
    if data.get("schema") != 4 or not str(data["version"]).startswith("8."):
        raise ValueError("Unsupported antiSMASH format; expected version 8, schema 4.")
    if checksum(genome["genbank_path"]) != genome["checksum"]:
        raise ValueError("GenBank no longer matches the stage-1 checksum.")
    with Path(genome["genbank_path"]).open(encoding="utf-8") as handle:
        originals = SeqIO.to_dict(SeqIO.parse(handle, "genbank"))
    evidence = {"regions": [], "gcf_matches": [], "issues": [], "provenance": {
        "genome": genome, "results_dir": str(results_dir), "result_json": str(result_path),
        "antismash_version": data["version"], "detector_settings": {}}}
    seen = set()
    try:
        for record in data["records"]:
            sequence_id = record.get("original_id") or record["id"]
            if sequence_id not in originals or sequence_id in seen:
                raise ValueError(f"Unexpected or duplicate result sequence: {sequence_id}")
            if record["seq"]["data"].upper() != str(originals[sequence_id].seq).upper():
                raise ValueError(f"Result sequence differs from the input: {sequence_id}")
            seen.add(sequence_id)
            _regions(record, sequence_id, genome["checksum"], data["version"], results_dir, evidence)
    except (KeyError, TypeError, IndexError) as exc:
        raise ValueError(f"Malformed antiSMASH structured results: {exc}") from exc
    if seen != set(originals):
        raise ValueError("antiSMASH results do not contain every input sequence.")
    return evidence


def _regions(record, sequence_id, digest, version, results_dir, evidence):
    detector = record["modules"].get("antismash.detection.hmm_detection", {})
    settings = {"strictness": detector.get("strictness"),
                "enabled_types": sorted(detector.get("enabled_types", [])),
                "multipliers": detector.get("rule_results", {}).get("multipliers")}
    evidence["provenance"]["detector_settings"][sequence_id] = settings
    qualifiers = {}
    for feature in record["features"]:
        if feature["type"] == "region":
            details = feature["qualifiers"]
            qualifiers[int(details["region_number"][0])] = details
    sequence = record["seq"]["data"]
    for number, area in enumerate(record["areas"], start=1):
        start, end = area["start"], area["end"]
        if not (0 <= start < len(sequence) and 0 < end <= len(sequence)) or start == end:
            raise ValueError(f"Invalid region coordinates in {sequence_id}.")
        segments = [[start, end]] if start < end else [[start, len(sequence)], [0, end]]
        details = qualifiers.get(number, {})
        identity = [digest, sequence_id, segments, version, settings, details.get("rules", [])]
        region_id = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        files = list(results_dir.rglob(f"{record['id']}.region{number:03d}.gbk"))
        if len(files) != 1:
            raise ValueError(f"Expected one region GenBank for {sequence_id}, region {number}.")
        with files[0].open(encoding="utf-8") as handle:
            bgc = SeqIO.read(handle, "genbank")
        expected = "".join(sequence[left:right] for left, right in segments)
        if str(bgc.seq).upper() != expected.upper():
            raise ValueError(f"BGC sequence differs from its region: {files[0]}")
        if not any(f.type == "region" for f in bgc.features):
            raise ValueError(f"Missing antiSMASH region annotation: {files[0]}")
        edge = details.get("contig_edge", [None])[0]
        evidence["regions"].append({
            "region_id": region_id, "sequence_id": sequence_id, "region_number": number,
            "start": start, "end": end, "segments": segments, "classes": area["products"],
            "crosses_origin": start > end, "contig_edge": None if edge is None else edge == "True",
            "bgc_path": str(files[0]), "bgc_checksum": checksum(files[0]),
        })
