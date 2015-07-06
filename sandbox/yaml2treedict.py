# -*- coding: utf-8 -*-

import yaml
from treedict import TreeDict


def yaml2treedict():

    with open("yaml2treedict.yaml", 'r') as f:
        content = yaml.load(f)
    config = TreeDict.fromdict(content, expand_nested=True)
    print config.makeReport()


if __name__ == "__main__":
    yaml2treedict()
