#!/usr/bin/env python
"""Script for creating and evil scan file from a series of images
"""
import os
import argparse

from pymapper import Scan


parser = argparse.ArgumentParser(description='Convert, image series to (evilish) scan file.')
parser.add_argument("plate", type=int,
    help="plate id")
parser.add_argument("imageDir", type=str,
    help="specify path to image files. laser1 and laser3 subdirectories are expected to exist. scan file will be placed here")
parser.add_argument("motorspeed", type=float,
    help="speed of motor in mm/s")
parser.add_argument("fps", type=float,
    help="frame rate in frames per second of camera")

if __name__ == "__main__":
    args = parser.parse_args()
    if not os.path.exists(args.imageDir):
        raise RuntimeError("couldn't locate image directory: %s"%args.imageDir)
    laser1dir = os.path.join(args.imageDir, "laser1")
    laser3dir = os.path.join(args.imageDir, "laser3")
    for laserDir in [laser1dir, laser3dir]:
        if not os.path.exists(laserDir):
            raise RuntimeError("couldn't locate image directory: %s"%laserDir)
    # create an outfilename
    fileName = "fiberScan-%i-py.par"%args.plate
    filePath = os.path.join(args.imageDir, fileName)
    scan = Scan(args.plate, filePath, motorSpeed=args.motorspeed, fps=args.fps)
    scan.batchProcess(laser1dir, motorID=1)
    scan.batchProcess(laser3dir, motorID=3)
    scan.writeOutputFile()
