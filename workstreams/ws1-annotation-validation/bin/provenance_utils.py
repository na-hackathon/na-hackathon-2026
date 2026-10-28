"""Helpers shared by ws1-provenance and ws1-manifest."""
import hashlib
import subprocess


def git(directory, *args):
    """Output of `git -C directory args`, or None when git is missing or the command fails."""
    try:
        result = subprocess.run(["git", "-C", str(directory), *args], capture_output=True, text=True)
    except OSError:
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip() or None


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()
