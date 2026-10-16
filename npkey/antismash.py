"""Stage 2: submit/resume antiSMASH or import its saved results."""

from html import unescape
import json
from pathlib import Path
import re
import time
from urllib.parse import quote, urljoin
from zipfile import ZipFile
import requests

from .antismash_results import parse_results
from .io_utils import checksum, get, write_json

BASE_URL = "https://antismash.secondarymetabolites.org"
SETTINGS = {"jobtype": "antismash8", "taxon": "bacteria", "genefinder": "none",
            "hmmdetection_strictness": "relaxed", "knownclusterblast": "false"}


def run_antismash(genome, output_dir, saved_results=None, *, poll_interval=30, timeout=3600):
    if poll_interval <= 0 or timeout < 0:
        raise ValueError("poll_interval must be positive and timeout non-negative.")
    run_dir = Path(output_dir).resolve() / "antismash" / genome["checksum"]
    run_dir.mkdir(parents=True, exist_ok=True)
    job = {}
    if saved_results is not None:
        results_dir = Path(saved_results).expanduser().resolve()
    else:
        checkpoint = run_dir / "job.json"
        if not checkpoint.exists():
            _submit(genome, checkpoint)
        job = json.loads(checkpoint.read_text(encoding="utf-8"))
        if not job.get("job_id"):
            raise RuntimeError(f"Submission uncertain; recover the job ID in {checkpoint}.")
        if job["parameters"] != SETTINGS:
            raise ValueError("Settings changed; import saved_results or use a new output directory.")
        results_dir = run_dir / "results"
        deadline = time.monotonic() + timeout
        while not job.get("results_downloaded") or not results_dir.is_dir():
            status = get(f"{BASE_URL}/api/v1.0/status/{quote(job['job_id'], safe='')}").json()
            job["status"] = status
            write_json(checkpoint, job)
            if status["state"] == "failed":
                raise RuntimeError(f"antiSMASH failed: {status.get('status')}")
            if status["state"] == "done":
                _download(urljoin(BASE_URL, status["result_url"]), run_dir)
                job["results_downloaded"] = True
                write_json(checkpoint, job)
                break
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError(f"Job {job['job_id']} still running; call again to resume.")
            time.sleep(min(poll_interval, remaining))
    evidence = parse_results(results_dir, genome)
    evidence["provenance"].update(job_id=job.get("job_id"), submission_parameters=job.get("parameters"))
    write_json(run_dir / "evidence.json", evidence)
    return evidence


def _submit(genome, checkpoint):
    path = Path(genome["genbank_path"])
    if checksum(path) != genome["checksum"]:
        raise ValueError("GenBank no longer matches the stage-1 checksum.")
    job = {"parameters": SETTINGS, "input_checksum": genome["checksum"], "job_id": None}
    with checkpoint.open("x", encoding="utf-8") as handle:
        json.dump(job, handle)
    try:
        with path.open("rb") as handle:
            response = requests.post(f"{BASE_URL}/api/v1.0/submit", data=SETTINGS,
                                    files={"seq": ("genome.gbk", handle, "application/octet-stream")},
                                    timeout=120, allow_redirects=False)
        response.raise_for_status()
        job["job_id"] = response.json()["id"]
        if not isinstance(job["job_id"], str) or not job["job_id"]:
            raise ValueError("Missing job ID")
    except (requests.RequestException, ValueError, KeyError) as exc:
        raise RuntimeError(f"Submission not confirmed; recover job ID in {checkpoint}.") from exc
    write_json(checkpoint, job)


def _download(url, run_dir):
    links = set(re.findall(r'href=["\']([^"\']+\.zip)["\']', get(url).text))
    if len(links) != 1:
        raise ValueError("Expected one antiSMASH ZIP download link.")
    archive_path = run_dir / "results.zip"
    archive_path.write_bytes(get(urljoin(url, unescape(links.pop()))).content)
    extracted = (run_dir / "results").resolve()
    with ZipFile(archive_path) as archive:
        for member in archive.infolist():
            if not (extracted / member.filename).resolve().is_relative_to(extracted):
                raise ValueError("Archive path outside its directory.")
        archive.extractall(extracted)
