

from loguru import logger
from pysiral.errorhandler import ErrorStatus


class DefaultLoggingClass(object):

    def __init__(self, cls_name):
        self.log = logger
        self.error = ErrorStatus(cls_name)