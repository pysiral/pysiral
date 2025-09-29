# -*- coding: utf-8 -*-

"""
This module defines default command line arguments for pysiral scripts.
"""

__author__ = "Stefan Hendricks <stefan.hendricks@awi.de>"

import argparse
from typing import List, Union, Callable, Any, Tuple, ClassVar, Optional
from dataclasses import dataclass, asdict, field


from pysiral import psrlcfg
from pysiral.core.flags import (
    Hemispheres, PysiralProcessingLevels, ProductProcessingLevels,
    DurationType, DataRecordType
)
from pysiral.scripts._argparse_types import (
    file_type, positive_int_type, dir_type, doi_type,
    pysiral_grid_id_type, month_number_type
)
from pysiral.scripts._argparse_actions import period_conversion, pysiral_settings_action


@dataclass(kw_only=True)
class ArgparseArgumentsArgs:
    name_or_flags: ClassVar[list[str]] = []
    action: Union[str, Callable, None] = field(default=None)
    nargs: Optional[str] = field(default=None)
    const: Union[str, int, None] = field(default=None)
    default: Any = field(default=None)
    type: Optional[Callable] = field(default=None)
    choices: List[Any] | None = None
    required: Optional[bool] = field(default=None)
    help: str = field(default=None)
    metavar: Optional[str] = field(default=None)
    dest: str = field(default=None)

    def get(self) -> Tuple[List[str], dict[str, Any]]:
        """
        Returns the input arguments and keyword arguments to `parser.add_argument()`

        :return: A tuple containing a list of argument names or flags and a
            ictionary of keyword arguments
        """
        # Convert the `args` property to a dictionary and remove empty values`
        args_dict = asdict(self)
        args_dict = {k: v for k, v in args_dict.items() if v is not None}

        # Remove the `name_or_flags` key from the dictionary and return it separately
        return self.name_or_flags, args_dict

    def as_positional(self, argument_name: str) -> "ArgparseArgumentsArgs":
        """
        Redefines the argument to be a positional argument

        :param argument_name: The argparse parser to add the argument to.
        """
        self.name_or_flags.clear()
        self.name_or_flags.append(argument_name)
        self.dest = None
        self.required = None
        return self


@dataclass(kw_only=True)
class PlatformID(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-p", "--platform"]
    action: str = "store"
    dest: str = "platform"
    type: Callable = str
    choices: List[Any] = field(default_factory=lambda: psrlcfg.platform_ids)
    metavar: str = "<platform_id>"
    help: str = """
    Radar altimeter platform id as defined in pysiral (see `pysiral info --platforms`). This option is 
    required only if the processor configuration file is applicable to multiple platforms 
    (e.g. sentinel3a, sentinel3b, etc.). If the processor configuration file is platform-specific, 
    this option has no effect.
    """


@dataclass(kw_only=True)
class ProcessingPeriod(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["processing_period"]
    action: Callable = period_conversion()
    metavar: str = "<processing period>"
    help: str = """
    Period definition for processing, given as a string in the format "YYYY-MM[-DD][:YYYY-MM[-DD]]".
    If only one date is given, it will be interpreted as a period (e.g., "2023-01" for January 2023 and
    ("2023-01-01" for one day). If two colon-separated dates are given, they will be interpreted 
    as a start and end date or month.
    """


@dataclass(kw_only=True)
class Hemisphere(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-H", "--hemisphere"]
    action: str = "store"
    dest: str = "hemisphere"
    default: Any = Hemispheres.GLOBAL
    metavar: str = "|".join(Hemispheres.get_choices())
    choices: List[Any] = field(default_factory=lambda: Hemispheres.get_choices())
    type: Callable = str
    help: str = """
    Target hemisphere for processing. Options are 'global', 'nh', or 'sh'. 
    """


@dataclass(kw_only=True)
class Duration(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-d", "--duration"]
    action: str = "store"
    dest: str = "duration"
    metavar: str = "<duration>"
    choices: List[Any] = field(default_factory=lambda: DurationType.get_choices())
    type: Callable = str
    help: str = """
    Duration type for splitting the processing period into segments. 
    Options are 'daily', 'isoweekly', or 'monthly' and will be passed
    `dateperiod.Dateperiod().get_segments()`. If left empty, the 
    duration type will be inferred from the processing period definition.
    """


@dataclass(kw_only=True)
class ProductProcessingLevel(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-p", "--processing-level"]
    dest: str = "processing_level"
    metavar: str = "<processing_level>"
    choices: List[Any] = field(default_factory=lambda: ProductProcessingLevels.get_choices())
    type: Callable = str
    help: str = """
    Target processing level. Setting this argument will overwrite any processing
    levels defined in the input data or processor definition.
    """+f" Valid processing levels are: [{', '.join(ProductProcessingLevels.get_choices())}]"


@dataclass(kw_only=True)
class DataRecord(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-D", "--data-record"]
    dest: str = "data_record"
    metavar: str = f"<data record>"
    choices: List[Any] = field(default_factory=lambda: DataRecordType.get_choices())
    type: Callable = str
    help: str = """
    Target data record. Setting this argument will overwrite any value
    defined in the input data or processor definition.
    """+f" Valid data records are: [{', '.join(DataRecordType.get_choices())}]"


@dataclass(kw_only=True)
class ExcludeMonths(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-E", "--exclude-months"]
    nargs: str = "+"
    type: Callable = month_number_type
    dest: str = "exclude_months"
    metavar: str = "<month number>"
    help: str = "List of months to be excluded from processing, given as integers (1-12)."


@dataclass(kw_only=True)
class InputVersion(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-v", "--input-version"]
    action: str = "store"
    dest: str = "input_version"
    default: Any = "default"
    metavar: str = "v{minor}p{major}|default"
    type: Callable = str
    help: str = """
    Input version name (e.g., v1p0, v2p1, or 'default'). 
    This is used to identify the version of the input data.
    """


@dataclass(kw_only=True)
class SourceDatasetID(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-s", "--source-dataset-id"]
    action: str = "store"
    dest: str = "source_dataset_id"
    metavar: str = "<source_dataset_id>"
    type: Callable = str
    help: str = """
    Identifier of the source dataset to be used for processing, summarizing the platform, version, 
    and timeliness information. The source dataset ID must be specified in the local machine definition file
    ({pysiral-cfg-location}/local_machine_def.yaml, e.g. `root.l1b_repository.<platform>.<source_dataset_id>`)
    """


@dataclass(kw_only=True)
class L1PFile(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["l1p_file"]
    type: Callable = file_type(suffix=".nc")
    help: str = """
    Path to Level-1 Pre-Processed (l1p) input file for the Level-2 processor.
    """


@dataclass(kw_only=True)
class L1PSettings(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["l1p_settings"]
    type: Callable = str
    action: Callable = pysiral_settings_action(
        target="proc",
        level=PysiralProcessingLevels.LEVEL1
    )
    metavar: str = "<id|filepath>"
    help: str = """
    Identifier or file path to the Level-1 Pre-Processor definition file.
    This file contains the settings for the Level-1 processor. The default location
    for these files is `{pysiral-cfg-location}/proc/l1/`. The identifier is the filename without 
    the `.yaml` extension. E.g.`cryosat2_pds_ipf1e_v1p2` will be resolved to
    `{pysiral-cfg-location}/proc/l1/cryosat2_pds_ipf1e_v1p2.yaml`.
    """


@dataclass(kw_only=True)
class L2Settings(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["l2_settings"]
    type: Callable = str
    action: Callable = pysiral_settings_action(
        target="proc",
        level=PysiralProcessingLevels.LEVEL2
    )
    metavar: str = "<l2 settings id|filepath>"
    help: str = """
    Identifier or file path to the Level-2 Processor definition file.
    This file contains the settings for the Level-2 processor. The default location
    for these files is `{pysiral-cfg-location}/proc/l2/`. The identifier is the filename without 
    the `.yaml` extension. E.g.`awi_cryosat2_nh_v2p6_rep` will be resolved to
    `{pysiral-cfg-location}/proc/l2/awi/v2p6/awi_cryosat2_nh_v2p6_rep.yaml`.
    """


@dataclass(kw_only=True)
class L3Settings(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["l3_settings"]
    action: Callable = pysiral_settings_action(
        target="proc",
        level=PysiralProcessingLevels.LEVEL3
    )
    type: Callable = str
    metavar: str = "<l3 settings id|filepath>"
    help: str = """
    Identifier or file path to the Level-3 Processor definition file.
    This file contains the settings for the Level-2 processor. The default location
    for these files is `{pysiral-cfg-location}/proc/l3/`. The identifier is the filename without 
    the `.yaml` extension. E.g.`l3c_cci_v4p0` will be resolved to
    `{pysiral-cfg-location}/proc/l3/cci/l3c_cci_v4p0.yaml`.
    """


@dataclass(kw_only=True)
class L2Outputs(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-o", "--l2-output"]
    nargs: str = "+"
    dest: str = "l2_outputs"
    action: Callable = pysiral_settings_action(target="output", level=ProductProcessingLevels.LEVEL2_INTERMEDIATE)
    metavar: str = "<l2 output id|filepath>"
    help: str = """
    Identifier or file path of one ore several Level-2 output definition files.
    Each file contains the output definition for the l2/l2i dataformat. The default location
    for these files is `{pysiral-cfg-location}/output/l2/`. The identifier is the filename without 
    the `.yaml` extension. E.g.`l2i_cci_v4p0` will be resolved to
    `{pysiral-cfg-location}/output/l2i/cci/l2i_cci_v4p0.yaml`.
    """


@dataclass(kw_only=True)
class L2POutputs(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-o", "--l2p-output"]
    action: Callable = pysiral_settings_action(
        target="output",
        level=ProductProcessingLevels.LEVEL2_PREPROCESSED
    )
    metavar: str = "<l2p output id|filepath>"
    help: str = """
    Identifier or file path of one Level-2 output definition file.
    Each file contains the output definition for the l2/l2i dataformat. The default location
    for these files is `{pysiral-cfg-location}/output/l2/`. The identifier is the filename without 
    the `.yaml` extension. E.g.`l2i_cci_v4p0` will be resolved to
    `{pysiral-cfg-location}/output/l2i/cci/l2i_cci_v4p0.yaml`.
    """


@dataclass(kw_only=True)
class L3Outputs(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-o", "--l3-output"]
    nargs: str = "+"
    type: Callable = str
    action: Callable = pysiral_settings_action(
        level=PysiralProcessingLevels.LEVEL3,
        target="output"
    )
    metavar: str = "<l3 output id|filepath>"
    help: str = """
    Identifier or file path of one or several Level-3 output definition files.
    Each file contains the output definition for the l3 dataformat. The default location
    for these files is `{pysiral-cfg-location}/output/l3/`. The identifier is the filename without 
    the `.yaml` extension. E.g.`l3c_cci_v4p0` will be resolved to
    `{pysiral-cfg-location}/output/l3/cci/l3c_cci_v4p0.yaml`.
    """


@dataclass(kw_only=True)
class MultiProcesssingNumCores(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-m", "--multiprocessing-num-cores"]
    dest: str = "multiprocessing_num_cores"
    type: Callable = positive_int_type
    metavar: str = "<num_cores>"
    default: Any = None
    help: str = """
    Set the number of CPU cores to be used for multiprocessing. If not set, the 
    default value (derived from `multiprocessing.cpu_count()`) will be used. 
    NOTE: In some managed environments, the default value is not reliable, which 
    may lead to performance issues. In this case, it is recommended to set this
    value manually to the known number of available CPU cores.
    """


@dataclass(kw_only=True)
class UseMultiProcesssing(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["--multiprocessing"]
    dest: str = "use_multiprocessing"
    action: Any = argparse.BooleanOptionalAction
    default: Any = True
    help: str = """
    Flag to allow disabling multiprocessing. If set, the processor will run in single-threaded mode.
    Default is True, meaning multiprocessing is enabled. 
    (also see option -m/--multiprocessing-num-cores)
    """


@dataclass(kw_only=True)
class ForceL2DefRecordType(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["--force-l2def-record-type"]
    dest: str = "force_l2def_record_type"
    default: Any = False
    help: str = """
    By default, the Level-2 processor will use the record_type defined in the l1p input files. 
    If this flag is set, the processor will use the record_type defined in the Level-2 definition file
    (l2def) instead. This is useful if the l1p input files are not consistent with the l2def record_type, 
    e.g. in case of metadata errors or if data from another timeliness is to be used for filling gaps. 
    This flag should only be set in special cases. 
    """


@dataclass(kw_only=True)
class L1PFiles(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["l1p_files"]
    nargs: str = "+"
    type: Callable = file_type("*.nc")
    metavar: str = "<l1p filepath>"
    help: str = """
    Target Level-1P (l1p) file(s) to be processed by the Level-2 processor.
    """


@dataclass(kw_only=True)
class L2iDirectory(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-i", "--l2i-product-dir", "--l2-product-dir"]
    nargs: Optional[str] = "+"
    dest: str = "l2_product_directory"
    type: Callable = dir_type(ends_with=["l2i", "l2"])
    metavar: str = "<l2 directory>"
    help: str = """
    Target Level-2i (l2i) product directory where the Level-2 output files will be written.
    The l2i files need to be organized in `yyyy/mm/` subdirectory structure. 
    """


@dataclass(kw_only=True)
class L3Directory(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["-O", "--l3-product-directory"]
    dest: str = "l3_product_directory"
    type: Callable = dir_type(
        must_exist=False,
        ends_with=[
            ProductProcessingLevels.LEVEL3_COLLATED,
            ProductProcessingLevels.LEVEL3_SUPERCOLLATED
        ]
    )
    metavar: str = "<l3 directory>"
    help: str = """
    Target Level-3 product directory where the Level-3 output files will be written. 
    The last sub-directory must be a valid Level-3 processing level code [l3c|l3s]. 
    If not given, the directory will be inferred from the l2 input directory. 
    """


@dataclass(kw_only=True)
class L3Grid(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["l3_grid_id"]
    type: Callable = pysiral_grid_id_type
    metavar: str = "<grid id>"
    help: str = """
    The grid definition id for the Level-3 processor. 
    """+f" Valid grid ids are: [{', '.join(psrlcfg.get_setting_ids('grid'))}]"


@dataclass(kw_only=True)
class DOI(ArgparseArgumentsArgs):
    name_or_flags: ClassVar[list[str]] = ["--doi"]
    type: Callable = doi_type
    metavar: str = "<product doi>"
    help: str = """
    The DOI (Digital Object Identifier) to be written in the global attributes of the product files.
    """


class DefaultCommandLineArguments(object):

    def __init__(self):

        self._args = {

            # Mission id
            "mission": {
                "action": 'store',
                "dest": 'mission_id',
                "choices": psrlcfg.platform_ids,
                "required": True,
                "help": "pysiral recognized mission id"
            },

            # platform (same as mission, but proper name)
            "platform": {
                "action": 'store',
                "dest": 'platform',
                "choices": psrlcfg.platform_ids,
                "required": True,
                "default": None,
                "help": "pysiral recognized platform id"
            },

            # Default date parameter
            "date": {
                "action": "store",
                "dest": "stop_date",
                "nargs": "+",
                "type": int,
                "required": False,
                "help": 'list as year and month and day (optional)'
            },

            # Default date parameter
            "hemisphere": {
                "action": "store",
                "dest": "hemisphere",
                "choices": ["global", "north", "south"],
                "default": "global",
                "required": False,
                "help": 'hemisphere flag for processing)'
            },

            # List of month to exclude from monthly iterations
            "exclude-month": {
                "action": "store",
                "dest": "exclude_month",
                "nargs": "+",
                "type": int,
                "required": False,
                "default": None,
                "help": 'list of months to be excluded from processing'
            },

            # Flag that indicates if previous versions shall be removed
            # before processing / plotting etc.
            "remove-old": {
                "action": "store_true",
                "dest": "remove_old",
                "default": False,
                "required": False,
                "help": 'remove all existing product in target directory'
            },

            # version tag of input data
            "input-version": {
                "action": "store",
                "dest": "input_version",
                "default": "default",
                "required": False,
                "help": 'input version name (see documentation)'
            },

            # same as input-version, but better worded
            "source-repo-id": {
                "action": "store",
                "dest": "source_repo_id",
                "default": None,
                "required": False,
                "help": 'specific tag in local_machine_def.yaml (root.l1b_repository.<platform>.<source_repo_od>'
            },

            # override any critical prompts for cronjobs etc
            "no-critical-prompt": {
                "action": "store_true",
                "dest": "no_critical_prompt",
                "default": False,
                "required": False,
                "help": 'set to skip any required command line inputs'
            },

            # preset for level-1b (l1bdata) fiels
            "l1b_files": {
                "action": "store",
                "dest": "l1b_files_preset",
                "default": None,
                "required": False,
                "help": 'Path to one or many l1bdata files (e.g.: path/*.nc)'
            },

            # fetch the level-1p file version
            "l1p-version": {
                "action": "store",
                "dest": "l1p-version",
                "default": None,
                "required": False,
                "help": 'file version of the l1p file'
            },

            # fetch the level-2 settings file
            "l1p-settings": {
                "action": "store",
                "dest": "l1p_settings",
                "default": None,
                "required": True,
                "help": 'id or path to Level-1P processor definition file file'
            },

            # fetch the level-2 settings file
            "l2-settings": {
                "action": "store",
                "dest": "l2_settings",
                "default": None,
                "required": True,
                "help": 'id or path to Level-2 settings file'
            },

            # fetch the level-2 settings file
            "l2-output": {
                "action": "store",
                "dest": "l2_output",
                "default": "l2i_default",
                "required": False,
                "help": 'l2 outputdef id'
            },

            # fetch the level-2 settings file
            "l2p-output": {
                "action": "store",
                "dest": "l2p_output",
                "default": "l2p_default",
                "required": False,
                "help": 'l2p outputdef id'
            },

            # set the run tag for the Level-2 Processor
            "run-tag": {
                "action": "store",
                "dest": "run_tag",
                "default": None,
                "required": False,
                "help": 'tag for the Level-2 output'
            },

            # no overwrite protection for level-2 outputs
            "no-overwrite-protection": {
                "action": "store_false",
                "dest": "overwrite_protection",
                "default": False,
                "required": False,
                "help": 'disable writing Level-2 output to unique directory'
            },

            # no overwrite protection for level-2 outputs
            "overwrite-protection": {
                "action": "store_true",
                "dest": "overwrite_protection",
                "default": False,
                "required": False,
                "help": 'enable writing Level-2 output to unique directory (default)'
            },

            "period": {
                "action": "store",
                "dest": "period",
                "default": "month",
                "required": False,
                "help": 'data period tag (default: month)'
            },

            "l2i-product-dir": {
                "action": "store",
                "dest": "l2i_product_dir",
                "nargs": "+",
                "default": None,
                "required": True,
                "help": "l2i input directory"
            },

            "l3-product-dir": {
                "action": "store",
                "dest": "l3_product_dir",
                "default": None,
                "required": False,
                "help": "l3 output directory"
            },

            "l3-settings": {
                "action": "store",
                "dest": "l3_settings",
                "default": "l3_default",
                "required": False,
                "help": "l3 settings definition id or filename"
            },

            "l3-griddef": {
                "action": "store",
                "dest": "l3_griddef",
                "default": None,
                "required": True,
                "help": "l3 grid definition id or filename"
            },

            "l3-output": {
                "action": "store",
                "dest": "l3_output",
                "default": "default",
                "required": True,
                "help": "l3 output id"
            },

            "doi": {
                "action": "store",
                "dest": "doi",
                "default": "None",
                "required": False,
                "type": str,
                "help": "doi number to be written in global attributes"
            },

            "data_record_type": {
                "action": "store",
                "dest": "data_record_type",
                "default": "None",
                "required": False,
                "type": str,
                "help": "type of data record [cdr, icdr]"
            },

            "force-l2def-record-type": {
                "action": "store_true",
                "dest": "force_l2def_record_type",
                "default": False,
                "required": False,
                "help": "overwrite l1p record type [cdr, icdr, nrt, ..] with metadata.record_type tag in l2def"
            },

            "mp-cpu-count": {
                "action": "store",
                "dest": "mp_cpu_count",
                "default": None,
                "type": int,
                "required": False,
                "help": "Number of CPU's to be used for multi-processing"
            },
        }

    def get_argparse_dict(self, name, destination, required):
        options = self._args[name]
        options["dest"] = destination
        options["required"] = required
        return options
