# -*- coding: utf-8 -*-

"""
=======
Purpose
=======

This module read the auto-generated version information toml files
and makes the importable.


"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import subprocess
try:
    import tomllib
except ImportError: # Python < 3.11
    import tomli as tomllib

from collections import defaultdict
from pathlib import Path
from typing import Union

_THIS_DIR = Path(__file__).parent.resolve()


def get_git_revision_hash(target="HEAD") -> Union[str, None]:
    try:
        return subprocess.check_output(['git', 'rev-parse', target], cwd=_THIS_DIR).decode('ascii').strip()
    except subprocess.SubprocessError:
        return None


def get_git_branch() -> Union[str, None]:
    try:
        return subprocess.check_output(['git', 'branch', '--show-current'], cwd=_THIS_DIR).decode('ascii').strip()
    except subprocess.CalledProcessError:
        return None

try:
    with open(_THIS_DIR / "software-version.toml", mode="rb") as fp:
        _software_version_data = tomllib.load(fp)
except FileNotFoundError:
    _software_version_data = defaultdict(None)


try:
    with open(_THIS_DIR / "git-version.toml", mode="rb") as fp:
        _git_version_data = tomllib.load(fp)
except FileNotFoundError:
    _software_version_data = defaultdict(None)

# Try to get the git information from local git, else default to
# git-version.toml
git_version_local = get_git_revision_hash()
git_branch_local = get_git_branch()

SOFTWARE_VERSION = _software_version_data.get("SOFTWARE_VERSION")
GIT_VERSION = git_version_local if git_version_local else _git_version_data.get("COMMIT")
GIT_BRANCH = git_branch_local if git_branch_local else _git_version_data.get("BRANCH")
GIT_ORIGIN = _git_version_data.get("ORIGIN")





