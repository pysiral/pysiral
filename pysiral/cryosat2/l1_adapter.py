
from pysiral.errorhandler import ErrorStatus
from pysiral.logging import DefaultLoggingClass


class ESAPDSBaselineD(DefaultLoggingClass):

    def __init__(self, cfg):
        cls_name = self.__class__.__name__
        super(ESAPDSBaselineD, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        self.cfg = cfg