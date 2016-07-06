#!/usr/bin/env python
import os
import argparse
import sys
import pickle

import numpy
from scipy.interpolate import UnivariateSpline

from pymapper.imgProcess import getImgTimestamps

class SlitHeadPosition(object):
    def __init__(self, imgTimestamps, start, speed):
        """
        !Callable class for determining slit head position from an image number or set of image numbers
        @param[in] imgTimestamps, list of timestamps corresponding to the image sequence
        @param[in] start: motor start position
        @param[in] speed: motor speed
        """
        imgNums = range(len(imgTimestamps))
        motorPositions = start + speed*imgTimestamps
        self.splineInterp = UnivariateSpline(
            x = imgNums,
            y = motorPositions,
            )

    def __call__(self, imgNum):
        """! solve for slithead position
        @param[in] imgNum float describing image number position (eg fractional image numbers allowed)
        @return motor position
        """
        return self.splineInterp(imgNum)


def extractValue(logLine):
    """Get the float value following the ":" in a line
    """
    return float(logLine.split(":")[-1].strip())

def getScanParams(logfile):
    """Parse the logfile to determine the scan params
    """
    speed = None
    start = None
    end = None
    with open(logfile, "r") as f:
        logLines = f.readlines()
    for line in logLines:
        if "motor start pos" in line:
            start = extractValue(line)
        elif "motor end pos" in line:
            end = extractValue(line)
        elif "motor scan speed" in line:
            speed = extractValue(line)
        if not None in [start, end, speed]:
            break
    if None in [start, end, speed]:
        raise RuntimeError("Could not extract start, end, and/or speed: (%s, %s, %s)"
            %(str(start), str(end), str(speed)))
    return start, end, speed


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Determine fiber positions on the slit."
        )
    parser.add_argument("scanDir", help="""Directory containing scan""")
    args = parser.parse_args()
    scanDir = os.path.abspath(args.scanDir)
    if not os.path.exists(scanDir):
        raise RuntimeError("Scan directory does not exit: %s"%scanDir)
    logfile = os.path.join(scanDir, "scan.log")
    if not os.path.exists(logfile):
        raise RuntimeError("Could not locate scan log: %s"%logfile)
    start, end, speed = getScanParams(logfile)
    # get the ordered detections
    detectionListFile = os.path.join(scanDir, "detectionList.pkl")
    if not os.path.exists(detectionListFile):
        raise RuntimeError("Could not locate detection list file: %s"%(detectionListFile))
    pkl = open(detectionListFile, "rb")
    detectedFiberList = pickle.load(pkl)
    pkl.close()
    if not len(detectedFiberList)==300:
        raise RuntimeError("Number of detections %i != 300"%len(detectedFiberList))
    imgTimestamps = getImgTimestamps(scanDir)
    slitHeadPosition = SlitHeadPosition(imgTimestamps, start, speed)
    fiberPositions = []
    for detection in detectedFiberList:
        # weight the detection frame by the counts
        imgNums = []
        counts = []
        for centroid in detection.centroidList:
            baseDir, imgFile = os.path.split(centroid["imageFile"])
            # img looks like img<int>.bmp, just get int part
            imgNums.append(int(imgFile[3:].split(".")[0])) # slice off img
            counts.append(centroid["counts"])
        weightedCenter = slitHeadPosition(numpy.average(imgNums, weights=counts))
        fiberPositions.append(weightedCenter)
    # write fiber positions to list
    fiberPosFile = os.path.join(scanDir, "fiberpos.dat")
    with open(fiberPosFile, "w") as f:
        f.writeLine("# Measured fiber positions for scan: %s"%scanDir)
        f.writeLine("# Fiber Number   Motor Position (mm)")
        for ind, fiberPos in enumerate(fiberPositions):
            f.writeLine("%i    %.6f"%(ind+1, fiberPos))

if __name__ == "__main__":
    main()
