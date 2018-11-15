
import os
from pysiral.errorhandler import ErrorStatus
from pysiral.logging import DefaultLoggingClass


class ESAPDSBaselineD(DefaultLoggingClass):

    def __init__(self, cfg, raise_on_error=False):

        cls_name = self.__class__.__name__
        super(ESAPDSBaselineD, self).__init__(cls_name)
        self.error = ErrorStatus(caller_id=cls_name)

        # Store arguments
        self.raise_on_error = raise_on_error
        self.cfg = cfg

        # Init main class variables
        self.nc = None

    def get_l1(self, filepath):
        """
        Main entry point to the CryoSat-2 Baseline-D Input Adapter
        :param filepath:
        :return:
        """

        # Input Validation
        if not os.path.isfile(filepath):
            msg = "Not a valid file: %s" % filepath
            self.log.warning(msg)
            self.error.add_error("invalid-filepath", msg)
            return self.empty

        # Parse the input file
        self._read_input_netcdf()




    @property
    def empty(self):
        return None