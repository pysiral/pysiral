# -*- coding: utf-8 -*-

"""
This submodule contains classes that mock retracker functionality.
"""

from loguru import logger

from pysiral.l2proc.procsteps import Level2ProcessorStep


class SGDRMultipleElevations(Level2ProcessorStep):
    """
    A processing step that computes elevation from a set of
    already retracked ranges and sets these as auxiliary
    variables

    NOTE: This step is only useful if the retracked ranges
          can be taken from a L2/SGDR file

          Also, if any corrections to the elevations are
          to be maded, the elevation parameter names
          have to be added as target_variables in the
          corresponding L2 processor steps
    """

    def __init__(self, *args, **kwargs):
        super(SGDRMultipleElevations, self).__init__(*args, **kwargs)

    def execute_procstep(self, l1b, l2):
        """
        Computes elevations from pre-defined retrackers
        :param l1b:
        :param l2:
        :return: error_status_flag
        """

        # Get a clean error status
        error_status = self.get_clean_error_status(l2.n_records)

        # Get options (and break if incorrect)
        classifier_name_fmt = self.cfg.options.get("classifier_name_fmt", None)
        output_name_fmt = self.cfg.options.get("output_name_fmt", None)
        if None in [classifier_name_fmt, output_name_fmt]:
            logger.error(f"- Invalid options to {self.__class__.__name__}, skipping")
            error_status[:] = True
            return error_status

        # Get initial elevations (altitude - range)
        # NOTE: It is assumed here that instrument and center of gravity correction
        #       have already been applied to the ranges in the l1 data
        for target_retracker in self.target_retrackers:

            # Get the retracked range from the l1 classifier data group
            classifer_var_name = classifier_name_fmt.format(target_retracker)
            retracker_range = l1b.classifier.get_parameter(classifer_var_name)

            # Compute elevation and add to l2
            elevation_value = l2.altitude[:] - retracker_range
            aux_id = self.auxid_fmt.format(target_retracker)
            aux_name = output_name_fmt.format(target_retracker)
            l2.set_auxiliary_parameter(aux_id, aux_name, elevation_value)

        # Return clean error status (for now)
        return error_status

    @property
    def target_retrackers(self):
        return self.cfg.options.get("predefined_retrackers", [])

    @property
    def l2_input_vars(self):
        return []

    @property
    def l2_output_vars(self):
        return [self.auxid_fmt.format(iter(self.target_retrackers))]

    @property
    def auxid_fmt(self):
        return "elev{}"

    @property
    def error_bit(self):
        return self.error_flag_bit_dict["other"]

