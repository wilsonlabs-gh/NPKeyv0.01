"""Stage 4: reference BGCs, documented compounds and explicit-CID InChIKeys."""

from datetime import datetime, timezone
import json
from pathlib import Path
import re
import tarfile
import time
import requests
from .io_utils import checksum, get, write_json


def add_compounds(evidence, metadata_archive, output_dir):
    model = evidence["provenance"]["gcf"]
    archive_path = Path(metadata_archive).expanduser().resolve()
    digest = checksum(archive_path)
    links = [(region["region_id"], accession) for region in evidence["regions"]
             if region["gcf_id"] is not None for accession in model["families"][str(region["gcf_id"])]]
    records = {}
    with tarfile.open(archive_path) as archive:
        for accession in sorted({accession for _, accession in links}):
            member = f"mibig_json_{model['mibig_release']}/{accession}.json"
            records[accession] = json.load(archive.extractfile(member))
            if records[accession]["accession"] != accession:
                raise ValueError(f"MIBiG accession mismatch: {member}")
    cache = Path(output_dir).resolve() / "pubchem"
    cache.mkdir(parents=True, exist_ok=True)
    lookups = {}
    evidence.update(reference_bgcs={}, compounds={}, bgc_compounds=[], region_bgcs=links)
    for accession, record in records.items():
        evidence["reference_bgcs"][accession] = {key: record.get(key) for key in ("version", "status", "quality")}
        for index, compound in enumerate(record.get("compounds", [])):
            cids = sorted({value.split(":")[1] for value in compound.get("databaseIds", [])
                           if re.fullmatch(r"pubchem:[0-9]+", value)})
            for cid in cids:
                if cid not in lookups:
                    lookups[cid] = _inchikey(cid, cache)
            keys = {lookups[cid].get("inchikey") for cid in cids}
            key = keys.pop() if len(keys) == 1 else None
            identity = key or f"mibig:{digest}:{accession}:{index}"
            evidence["compounds"][identity] = key
            evidence["bgc_compounds"].append([accession, identity, index, compound])
    evidence["provenance"]["compounds"] = {"mibig_release": model["mibig_release"],
        "archive": str(archive_path), "sha256": digest, "pubchem": lookups}
    return evidence


def _inchikey(cid, cache):
    path = cache / f"{cid}.json"
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChIKey/JSON"
    try:
        if path.exists():
            data = json.loads(path.read_text(encoding="utf-8"))
        else:
            time.sleep(0.21)
            data = dict(get(url).json(), url=url, retrieved_at=datetime.now(timezone.utc).isoformat())
        properties = data["PropertyTable"]["Properties"]
        if len(properties) != 1 or properties[0]["CID"] != int(cid):
            raise ValueError("PubChem CID mismatch")
        key = properties[0]["InChIKey"]
        if not isinstance(key, str) or not re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", key):
            raise ValueError("Invalid PubChem InChIKey")
        if not path.exists():
            write_json(path, data)
        return {"inchikey": key, "snapshot": str(path), "sha256": checksum(path),
                "url": data["url"], "retrieved_at": data["retrieved_at"]}
    except (requests.RequestException, ValueError, KeyError, TypeError) as exc:
        return {"url": url, "error": str(exc)}
