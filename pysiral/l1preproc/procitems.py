# -*- coding: utf-8 -*-

"""

"""

from schema import Schema, And, Or, Use
from typing import Any, Dict, TypeVar

from pysiral import get_cls
from pysiral.l1bdata import Level1bData


class L1PProcItem(object):

    def __init__(self, **cfg: Dict) -> None:
        """
        Base class for Level 1 Pre-Processor items. Not to be called direclty

        :param cfg: Configuration dictionary
        """
        self.cfg = cfg

    def apply(self, l1: Level1bData):
        """
        This class
        :param l1:
        :return:
        """
        raise NotImplementedError(f"Do not call {self.__class__.__name__} directly")

    def __getattr__(self, item: str) -> Any:
        """
        Direct attribute access to the cfg dictionary

        :param item:
        :return:
        """
        if item in self.cfg:
            return self.cfg[item]
        else:
            raise ValueError(f"attribute {item} not found in class or config dictionary")


L1PPROCITEM_CLS_TYPE = TypeVar("L1PPROCITEM_CLS_TYPE", bound=L1PProcItem)


class L1PProcItemDef(object):
    """
    Class for validating and processing Level-1 Pre-Processor item definition
    """

    def __init__(self,
                 label: str,
                 stage: str,
                 module_name: str,
                 class_name: str,
                 options_dict: dict = None
                 ) -> None:
        """
        Small helper class to digest the definition dictionary for
        Level-1 pre-processor items. The input are validated upon
        passing to this class.

        :param label: A label describing the processor item
        :param stage: The stage of the Level-1 processor where the item shall be executed
        :param module_name: The module name of the class (must be in local namespace)
        :param class_name: The specific class name
        :param options_dict:
        """

        self.label = Schema(str).validate(label)
        self.stage = Schema(And(str, lambda x: x in self.valid_stages)).validate(stage)
        self.module_name = Schema(str).validate(module_name)
        self.class_name = Schema(str).class_name
        option_dict = options_dict if options_dict is not None else {}
        self.option_dict = Schema(dict).validate(option_dict)
        breakpoint()

    @property
    def valid_stages(self):
        return ["post_source_file", "post_ocean_segment_extraction", "post_merge"]


def get_l1_proc_item(module_name: str,
                     class_name: str,
                     **cfg: Dict
                     ) -> L1PPROCITEM_CLS_TYPE:
    """
    A function returning an initialized processor item class

    :param module_name:
    :param class_name:
    :param cfg:

    :return: The initialized processor item class

    :raises: ValueError
    """

    # Get the class
    cls = get_cls(module_name, class_name, relaxed=True)

    # Check 1: class must exist
    if cls is None:
        raise ValueError(f"Unknown class {module_name}.{class_name} - Terminating")

    # Check 2: Class must be of inheriting L1PProcItem
    instance_ = cls(**cfg)
    parent_class = instance_.__class__.__bases__[0].__name__
    if parent_class != "L1PProcItem":
        raise ValueError(f"Class {module_name}.{class_name} does not bases on L1PProcItem ({parent_class})")

    return instance_
