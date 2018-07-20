# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 20:52:32 2016

@author: Stefan
"""

from pysiral.l1bpreproc import L1bPreProc
from pysiral.l1bdata import L1bConstructor
from pysiral.sentinel3.iotools import Sentinel3FileList

from pysiral.classifier import CS2PulsePeakiness
from pysiral.waveform import TFMRALeadingEdgeWidth

import time


class Sentinel3PreProc(L1bPreProc):

    def __init__(self):

        super(Sentinel3PreProc, self).__init__(self.__class__.__name__)

    def _get_input_files_local_machine_def(self, time_range, version):
        """
        Gets the input data files from default data repository that is
        specified in `local_machine_def.yaml`
        """

        # Get the path defined in local_machine_def.yaml
        l1b_repository = self._pysiral_config.local_machine.l1b_repository
        s3_l1b_repository = l1b_repository[self.mission_id][version].sral

        # Get the list of files (SAR and SIN in chronological order)
        # for the case of this test only for one month of data
        s3_files = Sentinel3FileList()
        s3_files.folder = s3_l1b_repository
        s3_files.search(time_range)

        # Transfer list of l1b input files
        self._l1b_file_list = s3_files.sorted_list

    def _get_l1bdata_ocean_segments(self, filename):
        """ Returns the source Sentinel-3 l1b data as a list of
        pysiral.L1bdata objects. """

        # Read the header information firsts
        try:
            l1b = L1bConstructor(self._pysiral_config)
            l1b.mission = self.mission_id
            l1b.filename = filename
            l1b.get_header_info()
        except:
            return []

        # 1) Check if file has ocean data in the polar regions at all
        has_polar_ocean = self.region_is_arctic_or_antarctic_ocean(l1b)

        # 2) Check if desires region is matched
        if self._jobdef.hemisphere == "global" or l1b.info.hemisphere == "global":
            matches_region = True
        elif l1b.info.hemisphere == self._jobdef.hemisphere:
            matches_region = True
        else:
            matches_region = False

        if not has_polar_ocean or not matches_region:
            return []

        # Only now read the full data set
        try:
            t0 = time.time()
            l1b.construct()
            t1 = time.time()
            self.log.info("- Parsed source file in %.3g seconds" % (t1 - t0))
        except:
            self.log.warning(" - Error reading l1b file")
            return []

        # Sentinel-3a WAT STC/NTC files cover a half-orbit from pole to plot
        # (with gaps). This approach is adapted from pre-processing
        # the Envisat half-orbits
        l1b_list = []
        if l1b.info.timeliness != "NRT":
            segments = self.extract_polar_segments_from_halforbit(l1b)

        # the NRT timeliness product however comes in 10 minute granules
        # therefore the content will be either arctic or antarctic data
        else:
            segments = [l1b]

        # Threshold for splitting
        seconds_threshold = self._mdef.max_connected_files_timedelta_seconds

        for seg in segments:

            # Get polar ocean segments
            if seg is not None:

                # Trim non ocean margins
                l1b_roi = self.trim_non_ocean_data(seg)

                # Test for unlikely case of no ocean data
                if l1b_roi is None:
                    return []

                # Split orbits at larger non-ocean segments
                l1b_segments = self.split_at_large_non_ocean_segments(l1b_roi)

                # There are data gaps in the L2WAT files, therefore split
                # at timestamp discontinuities
                l1b_segments = self.split_at_time_discontinuities(
                        l1b_segments, seconds_threshold, trim_non_ocean=True)

                # Some minor gaps seem to remain in the data
                # -> bring to regular grid
                for l1b_segment in l1b_segments:
                    try:
                        l1b_segment.detect_and_fill_gaps()
                    except AttributeError:
                        stop

                # Add to l1b stack
                try:
                    l1b_list.extend(l1b_segments)
                except TypeError:
                    self.log.info("- no ocean data in file")

        # Compute waveform based classifiers
        # The reason why this is done here is that it is very inefficient to
        # compute pulse peakiness and leading edge width for the entire
        # half-orbits. Doing it here means that this is done only for the
        # final data subset

        for l1b in l1b_list:

            # Get waveform data
            wfm = l1b.waveform.power
            rng = l1b.waveform.range
            radar_mode = l1b.waveform.radar_mode
            is_ocean = l1b.surface_type.get_by_name("ocean").flag

            tick = time.clock()
            # Calculate the Peakiness (CryoSat-2 notation)

            pulse = CS2PulsePeakiness(wfm)
            l1b.classifier.add(pulse.peakiness, "peakiness")
            l1b.classifier.add(pulse.peakiness_r, "peakiness_r")
            l1b.classifier.add(pulse.peakiness_l, "peakiness_l")
            tock = time.clock()
            print "CS2PulsePeakiness completed in %.1f seconds" % (tock-tick)

            tick = time.clock()
            # Compute the leading edge width (requires TFMRA retracking)
            lew = TFMRALeadingEdgeWidth(rng, wfm, radar_mode, is_ocean)
            lew1 = lew.get_width_from_thresholds(0.05, 0.5)
            lew2 = lew.get_width_from_thresholds(0.5, 0.95)
            l1b.classifier.add(lew1, "leading_edge_width_first_half")
            l1b.classifier.add(lew2, "leading_edge_width_second_half")
            l1b.classifier.add(lew.fmi, "first_maximum_index")
            tock = time.clock()
            print "TFMRALeadingEdgeWidth completed in %.1f seconds" % (tock-tick)

        return l1b_list

    @property
    def mission_id(self):
        # mission_id can be sentinel3a, ...
        return self._jobdef.mission_id
