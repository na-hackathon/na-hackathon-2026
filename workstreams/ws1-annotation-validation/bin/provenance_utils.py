"""Helpers shared by ws1-provenance and ws1-manifest."""
import hashlib
import subprocess
from pathlib import Path


def git(directory, *args):
    """Output of `git -C directory args`, or None when git is missing or the command fails."""
    try:
        result = subprocess.run(["git", "-C", str(directory), *args], capture_output=True, text=True)
    except OSError:
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def checkout_state(directory):
    """State of the git checkout in `directory`: 'clean', 'modified' or 'unavailable'.

    `git -C directory` works on the enclosing repository when `directory` is not a checkout itself (a
    submodule that was not initialised is an empty directory of its parent), so the directory has to be the
    top level of its repository. 'modified' means changes to tracked files; untracked files (build products)
    are not counted.
    """
    top = git(directory, "rev-parse", "--show-toplevel")
    if top is None or Path(top).resolve() != Path(directory).resolve():
        return "unavailable"
    return "modified" if git(directory, "status", "--porcelain", "--untracked-files=no") else "clean"


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()
