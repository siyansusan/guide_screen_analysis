#!/usr/bin/env python

from __future__ import (absolute_import, division, print_function,
                        unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
                      chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import os
import ProgramName
from Interval import Interval
from Rex import Rex

rex = Rex()
from StreamSamReads import StreamSamReads
from SamReadGroup import SamReadGroup
from SamHspFactory import SamHspFactory
from SamHspClusterer import SamHspClusterer
from SamAnnotation import SamAnnotation
from Ahab_allpairs import AhabAllPairs


# =========================================================================
#                         class AllPairsAnalysis
# =========================================================================
class AllPairsAnalysis(AhabAllPairs):
    def __init__(self, configFile, outputDir):
        super().__init__(configFile)
        config = self.config
        dedup = config.lookupOrDie("DEDUPLICATE")
        dedup = True if dedup in ("True", "true", "T", "TRUE", "Yes", "yes") \
            else False
        self.dedup = dedup
        self.MIN_IDENTITY = float(config.lookupOrDie("MIN_IDENTITY"))
        self.MAX_REF_GAP = float(config.lookupOrDie("MAX_REF_GAP"))
        self.MIN_ALIGNABILITY = float(config.lookupOrDie("MIN_ALIGNABILITY"))
        self.ESTIMATE_FIRST_CUT_SITE = float(config.lookupOrDie("ESTIMATE_FIRST_CUT_SITE"))
        self.ESTIMATE_SECOND_CUT_SITE = float(config.lookupOrDie("ESTIMATE_SECOND_CUT_SITE"))
        CUT_SITE_file = config.lookupOrDie("FIRST_CUT_SITE")
        self.MAX_ANCHOR_DISTANCE = int(config.lookupOrDie("MAX_ANCHOR_DISTANCE"))
        self.MAX_READ_GAP = int(config.lookupOrDie("MAX_READ_GAP"))
        self.MIN_ALIGNED_PROPORTION = \
            float(config.lookupOrDie("MIN_ALIGNED_PROPORTION"))
        self.OUTPUT_DIR = outputDir
        self.prepareOutputFiles(self.OUTPUT_DIR)
        self.prepareDebuggingFiles(self.OUTPUT_DIR)
        self.CHROMS = set(("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                           "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                           "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                           "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))
        firstsite = []
        secondsite = []
        with open(CUT_SITE_file) as first_site:
            for line in first_site:
                (guide_name, guide_site, guide_sequence) = line.strip().split("\t")
                if guide_name.startswith("V_50"):
                    firstsite.append([guide_name, guide_site, guide_sequence])
                else:
                    secondsite.append([guide_name, guide_site, guide_sequence])
        self.FIRST_CUT_SITE_LIST = firstsite
        self.SECOND_CUT_SITE_LIST = secondsite

    def __del__(self):
        self.BIN_ON_TARGET_DELETION.close()
        self.BIN_ON_TARGET_IMPERFECT_DELETION.close()
        self.BIN_ON_TARGET_NO_EDIT.close()
        self.BIN_ON_TARGET_SHORT_INDELS.close()
        self.BIN_ON_TARGET_INVERSION.close()
        self.BIN_OFF_TARGET_EDIT.close()
        self.BIN_OFF_TARGET_NO_EDIT.close()
        self.BIN_MISC.close()
        self.BIN_FAILED_FILTER.close()

    def prepareOutputFiles(self, DIR):
        """
         Prepare and open files containing basic information about the read.
         """
        if not os.path.exists(DIR):
            os.system("mkdir -p " + DIR)
        if not rex.find("(.*)/$", DIR):
            DIR = DIR + "/"
        self.BIN_ON_TARGET_DELETION = open(DIR + "bin-on-target-deletion.txt", "wt")
        self.BIN_ON_TARGET_IMPERFECT_DELETION = \
            open(DIR + "bin-on-target-imperfect-deletion.txt", "wt")
        self.BIN_ON_TARGET_NO_EDIT = open(DIR + "bin-on-target-no-edit.txt", "wt")
        self.BIN_ON_TARGET_SHORT_INDELS = \
            open(DIR + "bin-on-target-short-indels", "wt")
        self.BIN_ON_TARGET_INVERSION = \
            open(DIR + "bin-on-target-inversion.txt", "wt")
        self.BIN_OFF_TARGET_EDIT = open(DIR + "bin-off-target-edit.txt", "wt")
        self.BIN_OFF_TARGET_NO_EDIT = open(DIR + "bin-off-target-no-edit.txt", "wt")
        self.BIN_MISC = open(DIR + "bin-misc.txt", "wt")
        self.BIN_FAILED_FILTER = open(DIR + "bin-failed-filter.txt", "wt")

    def prepareDebuggingFiles(self, DIR):
        """
        Prepare and open files with all information about the read and HSPs in the read.
        """
        if not os.path.exists(DIR):
            os.system("mkdir -p " + DIR)
        if not rex.find("(.*)/$", DIR):
            DIR = DIR + "/"
        self.DEBUG_ON_TARGET_DELETION = \
            open(DIR + "debug-on-target-deletion.txt", "wt")
        self.DEBUG_ON_TARGET_DELETION_STAT = \
            open(DIR + "debug-on-target-deletion-stat.txt", "wt")
        self.DEBUG_ON_TARGET_IMPERFECT_DELETION = \
            open(DIR + "debug-on-target-imperfect-deletion.txt", "wt")
        self.DEBUG_ON_TARGET_NO_EDIT = \
            open(DIR + "debug-on-target-no-edit.txt", "wt")
        self.DEBUG_ON_TARGET_SHORT_INDELS = \
            open(DIR + "debug-on-target-short-indels", "wt")
        self.DEBUG_ON_TARGET_INVERSION = \
            open(DIR + "debug-on-target-inversion.txt", "wt")
        self.DEBUG_OFF_TARGET_EDIT = open(DIR + "debug-off-target-edit.txt", "wt")
        self.DEBUG_OFF_TARGET_NO_EDIT = \
            open(DIR + "debug-off-target-no-edit.txt", "wt")
        self.DEBUG_MISC = open(DIR + "debug-misc.txt", "wt")
        self.DEBUG_FAILED_FILTER = open(DIR + "debug-failed-filter.txt", "wt")

    def processCases(self, anno):
        """
        This method processes each read based on how many HSPs are in the read.
        """
        numHSPs = anno.numHSPs()
        if numHSPs == 1:
            self.process1HSP(anno)
        elif numHSPs == 2:
            self.process2HSPs(anno)
        elif numHSPs == 3:
            self.process3HSPs(anno)
        else:
            self.processManyHSPs(anno)

    def process1HSP(self, anno):
        """
        Precondition: 1 HSP only.
        """
        hsp = anno.getHSPs()[0]
        interval = hsp.getRefInterval()

        if (interval.contains(self.ESTIMATE_FIRST_CUT_SITE) or
                interval.contains(self.ESTIMATE_SECOND_CUT_SITE)):
            if hsp.containsIndels():
                self.bin(anno, self.BIN_ON_TARGET_SHORT_INDELS)
                self.dump(anno, self.DEBUG_ON_TARGET_SHORT_INDELS)
            else:
                self.bin(anno, self.BIN_ON_TARGET_NO_EDIT)
                self.dump(anno, self.DEBUG_ON_TARGET_NO_EDIT)
            return
        if hsp.containsIndels():
            self.bin(anno, self.BIN_OFF_TARGET_EDIT)
            self.dump(anno, self.DEBUG_OFF_TARGET_EDIT)
        else:
            self.bin(anno, self.BIN_OFF_TARGET_NO_EDIT)
            self.dump(anno, self.DEBUG_OFF_TARGET_NO_EDIT)

    def process2HSPs(self, anno):
        """
        Precondition: 2 HSPs only.
        """
        if anno.allSameStrand():
            self.process2HSPsSameStrand(anno)

    def checkRefGapDeletion(self, anno):
        gaps = anno.getRefGaps()
        if len(gaps) != 1:
            return False
        gap = gaps[0]
        if gap.getLength() > self.MAX_REF_GAP:
            self.bin(anno, self.BIN_FAILED_FILTER)
            self.dump(anno, self.DEBUG_FAILED_FILTER)
            return False
        deletion_region = Interval(self.ESTIMATE_FIRST_CUT_SITE, self.ESTIMATE_SECOND_CUT_SITE)
        if not gap.overlaps(deletion_region):
            self.bin(anno, self.BIN_OFF_TARGET_EDIT)
            self.dump(anno, self.DEBUG_OFF_TARGET_EDIT)
            return False
        return True

    def readGapSmallerThan(self, anno, MAX):
        gaps = anno.getReadGapLengths()
        if len(gaps) == 0:
            return True
        return gaps[0] < MAX

    def findClosestCutSites(self, anno):
        """
        Find guides that are closest to the cutsite
        """
        HSPs = anno.getHSPs()
        hsp1 = HSPs[0]
        hsp2 = HSPs[1]
        FIRST_CUT_SITE = self.FIRST_CUT_SITE_LIST
        SECOND_CUT_SITE = self.SECOND_CUT_SITE_LIST
        interval1 = hsp1.getRefInterval()
        interval2 = hsp2.getRefInterval()
        final_first_cut = None
        final_second_cut = None

        first_cut = 100000
        for site in FIRST_CUT_SITE:
            dis = int(site[1]) - interval1.getEnd()
            if 0 <= dis < first_cut:
                first_cut = dis
                final_first_cut = site
        second_cut = 100000
        for site in SECOND_CUT_SITE:
            dis = interval2.getBegin() - int(site[1])
            if 0 <= dis < second_cut:
                second_cut = dis
                final_second_cut = site

        return (final_first_cut, final_second_cut, str(first_cut), str(second_cut))

    def process2HSPsSameStrand(self, anno):
        """
        Precondition: 2 HSPs only and both HSPs are on the same strand.
        """
        HSPs = anno.getHSPs()
        hsp1 = HSPs[0]
        hsp2 = HSPs[1]
        interval1 = hsp1.getRefInterval();
        interval2 = hsp2.getRefInterval()
        if not self.checkRefGapDeletion(anno):
            return

        # Check some more filters
        if (not self.readGapSmallerThan(anno, self.MAX_READ_GAP) or
                anno.alignedProportion() < self.MIN_ALIGNED_PROPORTION):
            self.bin(anno, self.BIN_FAILED_FILTER)
            self.dump(anno, self.DEBUG_FAILED_FILTER)
            return

        (FIRST_CUT, SECOND_CUT, dis_1, dis_2) = self.findClosestCutSites(anno)
        if FIRST_CUT is None or SECOND_CUT is None:
            self.bin(anno, self.BIN_ON_TARGET_IMPERFECT_DELETION)
            self.dump(anno, self.DEBUG_ON_TARGET_IMPERFECT_DELETION)
            return
        else:

            FIRST_CUT_SITE = int(FIRST_CUT[1])
            SECOND_CUT_SITE = int(SECOND_CUT[1])
            guide_name = str(FIRST_CUT[0]) + "_" + str(SECOND_CUT[0])
            guide1 = str(FIRST_CUT[0])
            guide2 = str(SECOND_CUT[0])

            MAX_ANCHOR_DISTANCE = self.MAX_ANCHOR_DISTANCE
            d1 = FIRST_CUT_SITE - interval1.getEnd()
            d2 = interval2.getBegin() - SECOND_CUT_SITE
            if max(abs(d1), abs(d2)) > MAX_ANCHOR_DISTANCE:
                self.bin(anno, self.BIN_ON_TARGET_IMPERFECT_DELETION)
                self.dump_guide(anno, self.DEBUG_ON_TARGET_IMPERFECT_DELETION, guide_name)
                return
            self.bin(anno, self.BIN_ON_TARGET_DELETION)
            self.dump_guide(anno, self.DEBUG_ON_TARGET_DELETION, guide_name)
            self.dump_stat(anno, self.DEBUG_ON_TARGET_DELETION_STAT, guide1, guide2, dis_1, dis_2, d1, d2)

    def process3HSPs(self, anno):
        """
        Precondition: 3 HSPs only.
        """
        self.bin(anno, self.BIN_MISC)
        self.dump(anno, self.DEBUG_MISC)

    def processManyHSPs(self, anno):
        """
        Precondition: 4 HSPs or more.
        """
        self.bin(anno, self.BIN_MISC)
        self.dump(anno, self.DEBUG_MISC)

    def filter(self, anno):
        """
        This function applies some filters related to strand, target chromosome,
        percent identity, alignability, concordance, etc.
        """
        if not anno.allRefsSame():
            self.bin(anno, self.BIN_FAILED_FILTER)
            self.dump(anno, self.DEBUG_FAILED_FILTER)
            return False

        if anno.lowestPercentIdentity() < ahab.MIN_IDENTITY:
            self.bin(anno, self.BIN_FAILED_FILTER)
            self.dump(anno, self.DEBUG_FAILED_FILTER)
            return False

        if anno.anyRefsOverlap():
            self.bin(anno, self.BIN_FAILED_FILTER)
            self.dump(anno, self.DEBUG_FAILED_FILTER)
            return False
        return True


# =========================================================================
# main()
# =========================================================================
if len(sys.argv) != 4:
    exit(ProgramName.get() + " <settings.config> <filename.sam> <output-dir>\n")
(configFilename, samFile, outputDir) = sys.argv[1:]

# Instantiate Ahab object
ahab = AllPairsAnalysis(configFilename, outputDir)

# Process SAM file
hspFactory = SamHspFactory()
stream = StreamSamReads(samFile)
readsSeen = 0
read_dict = {}
while True:
    readGroup = stream.nextGroup()
    if not len(readGroup):
        break
    readsSeen += 1
    firstReads = readGroup.getReads()

    # Cluster the HSPs and produce an Alignment object
    HSPs = hspFactory.makeHSPs(firstReads)

    # dedup
    if len(HSPs) > 1:
        first_read_position = 100000000
        second_read_position = 100000000
        for hsp in HSPs:
            if hsp.getRec().flag_firstOfPair():
                if int(hsp.getRec().getRefPos()) < first_read_position:
                    first_read_position = int(hsp.getRec().getRefPos())
            if hsp.getRec().flag_secondOfPair():
                if int(hsp.getRec().getRefPos()) < second_read_position:
                    second_read_position = int(hsp.getRec().getRefPos())
        position = first_read_position, second_read_position
        if position in read_dict:
            continue
        else:
            read_dict[position] = HSPs[0].getReadID()

    HSPs = SamHspClusterer.cluster(HSPs)
    if len(HSPs) == 0: continue

    anno = SamAnnotation(HSPs)

    # Filter based on alignment quality and target chromosome, etc.
    if (not ahab.filter(anno)): continue

    # Address cases of 1 HSP, 2 HSPs, 3 HSPs, and >3 HSPs
    ahab.processCases(anno)
    del anno

print(ahab.readsBinned, "reads binned, out of ", readsSeen)
del ahab  # Call destructor to clean up
