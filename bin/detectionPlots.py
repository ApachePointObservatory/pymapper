"""recursive diving script for generating pyguide centroid plots.
"""
from __future__ import division, absolute_import

import traceback
import os

from pymapper.camera import unpickleCentroids, sortDectections

MINCOUNTS = 0
MINSEP = 3.5 # min between fibers separation in pixels

def findRawImgDirs():
    rawImgDirs = []
    for aDir in os.walk("/home/lcomapper/scan"):
        if "rawImage" in aDir[0]:
            rawImgDirs.append(aDir[0])
    return rawImgDirs


if __name__ == "__main__":
    imgDirs = findRawImgDirs()
    for imgDir in imgDirs:
        print("on dir", imgDir)
        try:
            sortDetections(unpickleCentroids(imgDir), plot=True)
        except:
            traceback.print_exc()
