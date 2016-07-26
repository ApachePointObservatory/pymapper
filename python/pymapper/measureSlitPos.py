#!/usr/bin/env python
import os
import argparse
import sys
import pickle

import numpy
from scipy.interpolate import UnivariateSpline

from .camera import getImgTimestamps, frameNumFromName, getScanParams

def getMeasuredFiberPositions(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    fiberPosList = []
    for line in lines:
        if line.strip().startswith("#"):
            continue # comment
        fiberNum, fiberPos = line.strip().split()
        fiberPosList.append(float(fiberPos))
    return numpy.asarray(fiberPosList)

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
    scanParams = getScanParams(logfile)
    start = scanParams["start"]
    speed = scanParams["speed"]
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
            imgNums.append(frameNumFromName(imgFile))
            counts.append(centroid["counts"])
        weightedCenter = slitHeadPosition(numpy.average(imgNums, weights=counts))
        fiberPositions.append(weightedCenter)
    # write fiber positions to list
    fiberPosFile = os.path.join(scanDir, "fiberpos.dat")
    with open(fiberPosFile, "w") as f:
        f.write("# Measured fiber positions for scan: %s\n"%scanDir)
        f.write("# Fiber Number   Motor Position (mm)\n")
        for ind, fiberPos in enumerate(fiberPositions):
            f.write("%i    %.6f\n"%(ind+1, fiberPos))

if __name__ == "__main__":
    main()
