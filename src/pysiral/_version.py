# -*- coding: utf-8 -*-

"""
=======
Purpose
=======

This module read the auto-generated version information toml files
and makes the importable.


"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

try:
    import tomllib
except ImportError: # Python < 3.11
    import tomli as tomllib

from collections import defaultdict
from pathlib import Path

_THIS_DIR = Path(__file__).parent.resolve()

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

SOFTWARE_VERSION = _software_version_data.get("SOFTWARE_VERSION")
GIT_VERSION = _git_version_data.get("COMMIT")
GIT_BRANCH = _git_version_data.get("BRANCH")
GIT_ORIGIN = _git_version_data.get("ORIGIN")





