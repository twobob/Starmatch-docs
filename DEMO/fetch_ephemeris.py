"""Download the official JPL DE200 ephemeris assets.

This repository cannot ship the binary ``de200.eph`` directly because the
code-hosting platform rejects non-text uploads. Running this helper will place
text-mode copies of the binary ephemeris and its header alongside the rest of
our demo data so the visualiser can remain entirely client-side.
"""

from __future__ import annotations

import hashlib
import shutil
import sys
from pathlib import Path
from typing import Iterable, Tuple
from urllib.request import urlopen

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
ASSETS: Tuple[Tuple[str, str, str], ...] = (
    (
        "de200.eph",
        "https://ephe.scryr.io/jpl/de200.eph",
        "1fa0594ee15717c924eccc302066eb151cc01aa4eb53e8f45afda98313fb7084",
    ),
    (
        "header.200",
        "https://ephe.scryr.io/jpl/header.200",
        "a7d9335afd2daac45c41e3be788f0ba86fef6ef4479118b78df69d24ee9f676b",
    ),
)


class DownloadError(RuntimeError):
    """Raised when a download fails or does not match the expected checksum."""


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ensure_data_dir() -> None:
    DATA_DIR.mkdir(exist_ok=True)


def download_file(dest: Path, url: str) -> None:
    try:
        with urlopen(url) as response, dest.open("wb") as target:
            shutil.copyfileobj(response, target)
    except OSError as exc:  # includes urllib errors
        raise DownloadError(f"Failed to download {url}: {exc}") from exc


def verify_file(dest: Path, checksum: str) -> None:
    if sha256sum(dest) != checksum:
        raise DownloadError(
            f"Checksum mismatch for {dest.name}. Delete the file and retry."
        )


def fetch_assets(assets: Iterable[Tuple[str, str, str]]) -> None:
    ensure_data_dir()
    for filename, url, checksum in assets:
        dest = DATA_DIR / filename
        if dest.exists():
            print(f"✓ {filename} already present")
            try:
                verify_file(dest, checksum)
            except DownloadError as exc:
                raise DownloadError(
                    f"Existing file {filename} does not match the expected checksum"
                ) from exc
            continue

        print(f"↓ Downloading {filename} …")
        download_file(dest, url)
        verify_file(dest, checksum)
        print(f"  Saved {dest.relative_to(DATA_DIR.parent)}")


def main() -> int:
    try:
        fetch_assets(ASSETS)
    except DownloadError as exc:
        print(exc, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
