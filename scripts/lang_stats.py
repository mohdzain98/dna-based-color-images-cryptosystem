"""Estimate GitHub's language bar (Python vs Jupyter Notebook) for this repo.

GitHub Linguist weighs languages by file size in bytes, so we do the same.
Usage: python scripts/lang_stats.py [repo_path]
"""

import subprocess
import sys
from pathlib import Path

LANGUAGES = {".py": "Python", ".ipynb": "Jupyter Notebook"}
SKIP_DIRS = {".git", ".venv", "venv", "__pycache__", ".pytest_cache", ".ruff_cache"}


def list_files(root: Path) -> list[Path]:
    """Return tracked and untracked (non-ignored) files via git, or walk the tree if git is unavailable."""
    try:
        out = subprocess.run(
            ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
            cwd=root,
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        return [root / line for line in out.splitlines()]
    except (subprocess.CalledProcessError, FileNotFoundError):
        return [
            p
            for p in root.rglob("*")
            if p.is_file() and not SKIP_DIRS.intersection(p.parts)
        ]


def count_bytes(files: list[Path]) -> dict[str, int]:
    """Sum file sizes per language, ignoring other file types."""
    totals = {name: 0 for name in LANGUAGES.values()}
    for path in files:
        language = LANGUAGES.get(path.suffix)
        if language and path.is_file():
            totals[language] += path.stat().st_size
    return totals


def main() -> None:
    """Print bytes and percentage per language."""
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    totals = count_bytes(list_files(root))
    grand_total = sum(totals.values()) or 1
    for language, size in sorted(totals.items(), key=lambda kv: -kv[1]):
        print(f"{language:<18}{size:>10,} bytes  {size / grand_total * 100:6.2f}%")


if __name__ == "__main__":
    main()
