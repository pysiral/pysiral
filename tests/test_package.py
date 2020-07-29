# -*- coding: utf-8 -*-
"""
@author: Stefan
"""

import unittest
import importlib
import pkgutil


class TestPackage(unittest.TestCase):

    def setUp(self):
        pass

    def testAllImports(self):
        import_submodules("pysiral")


def import_submodules(package, recursive=True):
    """ Import all submodules of a module, recursively, including subpackages

    :param package: package (name or actual module)
    :type package: str | module
    :rtype: dict[str, types.ModuleType]
    """
    if isinstance(package, str):
        package = importlib.import_module(package)
    results = {}
    for loader, name, is_pkg in pkgutil.walk_packages(package.__path__):
        full_name = package.__name__ + '.' + name
        # Make an exception for the cythonized part of pysiral
        if "bnfunc" in full_name:
            continue
        results[full_name] = importlib.import_module(full_name)
        if recursive and is_pkg:
            results.update(import_submodules(full_name))
    return results


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestPackage)
    unittest.TextTestRunner(verbosity=2).run(suite)
