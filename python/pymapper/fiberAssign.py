"""Routines for identifying, assigning fibers. Determining which (if any are missing)
"""
from __future__ import division, absolute_import

import os
import itertools

import numpy

import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt

def frameNumFromName(imgName, imgBase="img", imgExt="bmp"):
    return int(imgName.split(imgBase)[-1].split(".%s"%imgExt)[0])

if __name__ == "__main__":
    import pickle
    imgDir = "/home/lcomapper/Desktop/mappertestingMarch/run6"

    pkl_file = open(os.path.join(imgDir, "pickledDetections.pkl"), "rb")

    detectedFiberList = pickle.load(pkl_file)
    pkl_file.close()

    firstFrameNum = frameNumFromName(detectedFiberList[0]["imageFrames"][0])
    lastFrameNum = frameNumFromName(detectedFiberList[-1]["imageFrames"][-1])
    nDetectedFibers = len(detectedFiberList)

    # create visualization
    imgArray = numpy.zeros((nDetectedFibers, lastFrameNum-firstFrameNum+1))
    fig = plt.figure(figsize=(20,5))
    for fiberNum, detectedFiber in enumerate(detectedFiberList):
        for counts, filename in itertools.izip(detectedFiber["counts"], detectedFiber["imageFrames"]):
            frameNumber = frameNumFromName(filename)
            imgInd = frameNumber - firstFrameNum
            imgArray[fiberNum, imgInd] = counts / float(numpy.max(detectedFiber["counts"]))
            # plt.plot(fiberNum, imgInd, ".k")

    plt.imshow(imgArray, interpolation="none")
    nfn = os.path.join(imgDir, "fiberMap.png")
    plt.show()
    #fig.savefig(nfn)
    #plt.close(fig)    # close the figure
