
# TODO: Evaluate use and move to core module

from pysiral.core.errorhandler import ErrorStatus


class DefaultLoggingClass(object):
    """
    Template for default pysiral class with logging/error handling capabilities
    """

    def __init__(self, cls_name=None):
        """
        Init the class with a (loguru) logger and an ErrorStatus error handler
        :param cls_name:
        """

        self.error = ErrorStatus(cls_name)
