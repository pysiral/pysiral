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
        return subprocess.check_output(
            ['git', 'rev-parse', target],
            cwd=_THIS_DIR,
            stderr=subprocess.PIPE
        ).decode('ascii').strip()
    except subprocess.SubprocessError:
        return None


def get_git_branch() -> Union[str, None]:
    try:
        return subprocess.check_output(
            ['git', 'branch', '--show-current'],
            cwd=_THIS_DIR,
            stderr=subprocess.PIPE
        ).decode('ascii').strip()
    except subprocess.CalledProcessError:
        return None


def get_git_origin() -> Union[str, None]:
    try:
        return subprocess.check_output(
            ['git', 'remote', 'get-url', 'origin'],
            cwd=_THIS_DIR,
            stderr=subprocess.PIPE
        ).decode('ascii').strip()
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
GIT_VERSION_LOCAL = get_git_revision_hash()
GIT_BRANCH_LOCAL = get_git_branch()
GIT_ORIGIN_LOCAL = get_git_origin()

SOFTWARE_VERSION = _software_version_data.get("SOFTWARE_VERSION")
GIT_VERSION = GIT_VERSION_LOCAL if GIT_VERSION_LOCAL else _git_version_data.get("COMMIT")
GIT_BRANCH = GIT_BRANCH_LOCAL if GIT_BRANCH_LOCAL else _git_version_data.get("BRANCH")
GIT_ORIGIN = GIT_ORIGIN_LOCAL if GIT_ORIGIN_LOCAL else _git_version_data.get("ORIGIN")
