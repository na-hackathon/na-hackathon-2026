"""envs/ws1.lock.yml (exact version of every package) must agree with the direct dependencies of
envs/ws1.yml.

Run with:  python -m pytest tests/test_env_lock.py
"""
from __future__ import annotations

import re
from pathlib import Path

ENVS = Path(__file__).resolve().parents[1] / "envs"
TEST_TOOLS = {"pytest", "pydantic", "iniconfig", "pluggy"}      # installed next to the env to run tests


def packages(path: Path) -> dict[str, str]:
    """{normalised name: version or version pattern} of the conda and pip entries of an env file."""
    found = {}
    for line in path.read_text().splitlines():
        entry = re.match(r"\s+- (?:[\w.-]+::)?([A-Za-z0-9_.-]+?)(?:==|=)([^\s#]+)", line)
        if entry:
            found[entry.group(1).lower().replace("_", "-")] = entry.group(2)
    return found


def matches(pattern: str, version: str) -> bool:
    base = pattern.rstrip("*").rstrip(".")
    return version == base or version.startswith(base + ".")


def test_lock_contains_every_direct_dependency_at_the_pinned_minor_version():
    spec, lock = packages(ENVS / "ws1.yml"), packages(ENVS / "ws1.lock.yml")
    assert spec, "no pinned dependency read from envs/ws1.yml"
    for name, pattern in spec.items():
        assert name in lock, f"{name} missing from envs/ws1.lock.yml"
        assert matches(pattern, lock[name]), f"{name}: lock has {lock[name]}, envs/ws1.yml pins {pattern}"


def test_lock_pins_the_transitive_dependencies_and_not_the_test_tools():
    lock = packages(ENVS / "ws1.lock.yml")
    assert len(lock) > 2 * len(packages(ENVS / "ws1.yml"))          # more packages than the direct ones
    assert not TEST_TOOLS & lock.keys()
