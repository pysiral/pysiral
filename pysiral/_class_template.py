

from loguru import logger
from pysiral.errorhandler import ErrorStatus


class DefaultLoggingClass(object):
    """
    Template for default pysiral class with logging/error handling capabilities
    """

    def __init__(self, cls_name=None):
        """
        Init the class with a (loguru) logger and an ErrorStatus error handler
        :param cls_name:
        """
        self.log = logger
        self.error = ErrorStatus(cls_name)
