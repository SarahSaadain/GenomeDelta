import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "linux" / "scripts"
MAIN_SH = REPO_ROOT / "linux" / "main.sh"


def _run_python_script(name, args, cwd=None):
    """Run one of the linux/scripts/*.py CLI scripts as a subprocess."""
    script = SCRIPTS_DIR / name
    return subprocess.run(
        [sys.executable, str(script), *[str(a) for a in args]],
        cwd=cwd,
        capture_output=True,
        text=True,
    )


def _run_bash_script(script_path, args, cwd=None):
    return subprocess.run(
        ["bash", str(script_path), *[str(a) for a in args]],
        cwd=cwd,
        capture_output=True,
        text=True,
    )


@pytest.fixture
def run_script():
    """Run a script in linux/scripts/ by filename, e.g. run_script('MSA2consensus.py', [...])."""
    return _run_python_script


@pytest.fixture
def run_bash():
    """Run an arbitrary bash script, e.g. run_bash(MAIN_SH, [...])."""
    return _run_bash_script


@pytest.fixture
def main_sh():
    return MAIN_SH
