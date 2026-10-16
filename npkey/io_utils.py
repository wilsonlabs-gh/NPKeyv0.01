"""Shared file integrity, atomic JSON, read retries and external-command I/O."""

from contextlib import closing
import hashlib
import json
from pathlib import Path
import sqlite3
import subprocess
import time
import requests


def checksum(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path, data):
    path = Path(path)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def get(url):
    for attempt in range(3):
        try:
            response = requests.get(url, timeout=120)
            response.raise_for_status()
            return response
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError,
                requests.exceptions.ChunkedEncodingError) as exc:
            if isinstance(exc, requests.HTTPError):
                if exc.response.status_code not in {408, 429, 500, 502, 503, 504}:
                    raise
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def rows(path, sql, parameters=()):
    with closing(sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True)) as db:
        db.row_factory = sqlite3.Row
        return [dict(row) for row in db.execute(sql, parameters)]


def execute(command, log):
    with Path(log).open("w", encoding="utf-8") as handle:
        subprocess.run([str(arg) for arg in command], check=True, stdin=subprocess.DEVNULL,
                       stdout=handle, stderr=subprocess.STDOUT)
