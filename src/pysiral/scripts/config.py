# -*- coding: utf-8 -*-

"""

"""

from pathlib import Path
from typing import Literal, Optional

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"



def config(
    action: Literal["set", "update", "show"] = "show",
    target_directory: Optional[Path] = None,
    copy_local_machine_def_template: bool = False
) -> None:
    pass


