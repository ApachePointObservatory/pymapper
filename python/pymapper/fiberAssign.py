"""Routines for identifying, assigning fibers. Determining which (if any are missing)
"""
from __future__ import division, absolute_import

import sdss

import os
import itertools

import numpy
from scipy import stats

import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt

# spacing between blocks (measured empirically), units are adjacent fiber spacing
# turns out to be better to just use the same spacing for every block...
# BlockSpaceMultiplier = [1.928, 1.868, 1.828, 1.877, 1.910, 1.805, 1.980, 1.917, 1.887]

DEBUG = False

class ModeledTrace(object):
    def __init__(self, fiberSpacing, blockSpacing):
        """Model for a slithead trace
        """
        self.nBlocks = 10
        self.fibersPerBlock = 30
        self.fiberSpacing = fiberSpacing
        # 9 spaces for 10 blocks, should we allow each block to float a little bit? for now assume they're fixed
        self.blockSpacing = blockSpacing
        self.constructModelTrace()

    def constructModelTrace(self):
        maxFramesRequired = int(numpy.round((self.nBlocks-1)) * self.blockSpacing + (self.fibersPerBlock-1)*self.nBlocks * self.fiberSpacing)
        # add in extra frames as a buffer on either side of the array we will populate
        # these will be zeros
        bufferFrames = maxFramesRequired // 10
        startFrame = bufferFrames // 2
        maxFramesRequired += bufferFrames
        self.modeledTrace = numpy.zeros(maxFramesRequired)
        blockBoundaries = []
        for block in range(self.nBlocks):
            blockBegInd = None
            for nFiber in range(self.fibersPerBlock):
                position = block * (self.blockSpacing + (self.fibersPerBlock - 1) * self.fiberSpacing) + nFiber * self.fiberSpacing
                # set the nearest frame to 1 and its neighbors to 0.5
                # choose nearest frame because fiberSpacing and blockSpacing are floats
                posInt = int(numpy.round(position)) + startFrame
                # print("posInt", posInt, posInt1)
                if not numpy.all(self.modeledTrace[posInt-2:posInt+3]==0):
                    raise RuntimeError("Overlapping modeled slit detections")
                # construct a top-hat mask-like thing spanning 5 frames for each modeled detection, only the center gets
                # 1, so we can easily search for it later
                self.modeledTrace[posInt-2] = 0.99
                self.modeledTrace[posInt+2] = 0.99
                self.modeledTrace[posInt-1] = 0.99
                self.modeledTrace[posInt+1] = 0.99
                self.modeledTrace[posInt] = 1
                if blockBegInd is None:
                    # record block starting point
                    blockBegInd = posInt - 2
            # record block ending point
            blockEndInd = posInt + 2
            blockBoundaries.append((blockBegInd, blockEndInd))
        self.blockBoundaries = blockBoundaries

    def findMissingFibers(self, normalizedMeasuredTrace, missingFiberNum):
        # measured trace must have detections normalized (highest count detection = 1)
        nFramesMeas = len(normalizedMeasuredTrace)
        nFramesModel = len(self.modeledTrace) # remember model is buffered with zeros on either side
        nShifts = nFramesModel - nFramesMeas
        if nShifts < 0:
            raise RuntimeError("Measured trace has detetections too far separated in time")
        mostCounts = 0
        for ii in range(nShifts):
            # a sliding window, find the best position with most overlap between model and measure
            counts = numpy.sum(self.modeledTrace[ii:ii+nFramesMeas] * normalizedMeasuredTrace)
            if counts > mostCounts:
                mostCounts = counts
                bestFitIndex = ii

        # next figure out which fibers we got!!!
        allMissingFibers = []
        for nBlock in range(self.nBlocks):
            # determine the blockBoundaries for the measured trace
            blockBegIndModel, blockEndIndModel = self.blockBoundaries[nBlock]
            blockBegIndMeas = blockBegIndModel - bestFitIndex
            blockEndIndMeas = blockEndIndModel - bestFitIndex

            if blockBegIndMeas < 0:
                blockBegIndMeas = 0
            # if upper index is beyond length of measuredArray, doesn't matter for slicing purposes
            # how many detections did we find in this block
            foundInMeasuredBlock = numpy.argwhere(normalizedMeasuredTrace[blockBegIndMeas:blockEndIndMeas] == 1)
            if len(foundInMeasuredBlock) < self.fibersPerBlock:
                # not all fibers were found for this block
                # determine which ones are missing
                foundInModelBlock = numpy.argwhere(self.modeledTrace[blockBegIndModel:blockEndIndModel] == 1)
                assert len(foundInModelBlock) == self.fibersPerBlock, "must have 30 modeled detections. BUG!!!"
                matchedFiberPos = []
                for fiberPos in foundInMeasuredBlock:
                    dist = numpy.abs(fiberPos - foundInModelBlock)
                    # print("min dist", numpy.min(dist)/self.fiberSpacing)
                    if numpy.min(dist) > 0.5*self.fiberSpacing:
                        # amgiguous match?
                        raise RuntimeError("I think this could be an ambiguous match. Raising for paranoia for now.")
                    matchedFiberPos.append(numpy.argmin(dist))

                # pretty sure it would be impossible to duplicate matches, but just in case...
                if len(matchedFiberPos) != len(set(matchedFiberPos)):
                    raise RuntimeError("Detection matched two fibers in model")
                # lastly discover which fibers are missing from this block
                missingFibers = list(set(range(self.fibersPerBlock)) - set(matchedFiberPos))
                missingFibers = [nBlock*30 + fiberNum for fiberNum in missingFibers]
                # print("missing fibers in block!!!", nBlock, missingFibers)
                allMissingFibers.extend(missingFibers)
            # print("found ", len(foundInMeasuredBlock), "in block", nBlock+1)
        assert len(allMissingFibers) == missingFiberNum, "Didn't find same amount of missing fibers!!!!"

        # next allow individual blocks (groups of 30 fibers) to shift independently
        # determine the shift to apply to each block


        # we've found the best fit trace now identify the missing ones
        if DEBUG == True:
            plt.figure(figsize=(20,5))
            plt.plot(range(nFramesMeas), normalizedMeasuredTrace)
            plt.plot(numpy.arange(nFramesModel)-bestFitIndex, self.modeledTrace)
            plt.show()
        # print("bestFitIndex", bestFitIndex, "counts", mostCounts)

        return allMissingFibers



class FocalSurfaceSolver(object):
    def __init__(self, detectedFiberList, plPlugMap):
        pass

class SlitheadSolver(object):
    def __init__(self, detectedFiberList):
        """! detectedFiberList
        """
        self.detectedFiberList = detectedFiberList
        self.totalFiberNum = 300
        self.detectedFiberNum = len(detectedFiberList)
        self.missingFiberNum = self.totalFiberNum - self.detectedFiberNum
        self.firstFrameNum = frameNumFromName(detectedFiberList[0]["imageFrames"][0])
        self.lastFrameNum = frameNumFromName(detectedFiberList[-1]["imageFrames"][-1])
        self.totalFrames = self.lastFrameNum - self.firstFrameNum + 1
        # for each detection determine the center frame
        # weight each frame by the detected counts
        detectionCenters = []
        for ind, detectedFiber in enumerate(detectedFiberList):
            frameNumbers = [frameNumFromName(frameName) for frameName in detectedFiber["imageFrames"]]
            counts = detectedFiber["counts"]
            detectionCenters.append(numpy.average(frameNumbers, weights=counts))
        diffDetectionCenters = numpy.diff(detectionCenters)
        roughFiberSpacing = numpy.median(diffDetectionCenters)
        # throw out outliers and recompute.  This probably doesn't do much
        self.fiberSpacing = numpy.median(diffDetectionCenters[numpy.nonzero(diffDetectionCenters < roughFiberSpacing*1.3)])
        # for determining individual block spacing...must have all detections in all fibers...
        # print("block space 1", (detectionCenters[30]-detectionCenters[29])/self.fiberSpacing)
        # print("block space 2", (detectionCenters[60]-detectionCenters[59])/self.fiberSpacing)
        # print("block space 3", (detectionCenters[90]-detectionCenters[89])/self.fiberSpacing)
        # print("block space 4", (detectionCenters[120]-detectionCenters[119])/self.fiberSpacing)
        # print("block space 5", (detectionCenters[150]-detectionCenters[149])/self.fiberSpacing)
        # print("block space 6", (detectionCenters[180]-detectionCenters[179])/self.fiberSpacing)
        # print("block space 7", (detectionCenters[210]-detectionCenters[209])/self.fiberSpacing)
        # print("block space 8", (detectionCenters[240]-detectionCenters[239])/self.fiberSpacing)
        # print("block space 9", (detectionCenters[270]-detectionCenters[269])/self.fiberSpacing)
        self.blockSpacing = numpy.median(diffDetectionCenters[numpy.nonzero(diffDetectionCenters > self.fiberSpacing*1.3)])
        self.normalizedMeasuredTrace = self.getNormalizedMeasuredTrace()
        self.modeledTrace = ModeledTrace(self.fiberSpacing, self.blockSpacing)
        # get missing fibers
        self.missingFiberNumbers = self.modeledTrace.findMissingFibers(self.normalizedMeasuredTrace, self.missingFiberNum)
        # measuredTrace = self.getMeasuredTrace()
        # self.plotTrace(measuredTrace)
        # optimalTrace = self.getOptimalTrace()
        # self.plotTrace(optimalTrace)

    def getNormalizedMeasuredTrace(self):
        # construct from the detectedFiberList
        # Counts are normalized by max value (thus cannot exceed 1)
        # if multiple frames share a maximum value, force all except the first
        # to a lower value.  Allow only one to have the max value
        # this is for later searching
        measuredTrace = numpy.zeros(self.totalFrames)
        for detectedFiber in self.detectedFiberList:
            maxCounts = numpy.max(detectedFiber["counts"])
            foundMax = False
            for imgFrame, counts in itertools.izip(detectedFiber["imageFrames"], detectedFiber["counts"]):
                frameNumber = frameNumFromName(imgFrame)
                traceInd = frameNumber - self.firstFrameNum
                normalizedCounts = counts / maxCounts
                if normalizedCounts == 1:
                    if foundMax:
                        # only allow one detection to equal 1 (the first found)
                        # to ease in searching for detections later
                        normalizedCounts = 0.999
                    else:
                        foundMax = True
                assert measuredTrace[traceInd] == 0, "already counts found in this frame?!?!"
                measuredTrace[traceInd] = normalizedCounts
        return measuredTrace

def frameNumFromName(imgName, imgBase="img", imgExt="bmp"):
    return int(imgName.split(imgBase)[-1].split(".%s"%imgExt)[0])

if __name__ == "__main__":
    import pickle
    imgDir = "/Users/csayres/Desktop/mapperPyguide/run6noFlat"

    pkl_file = open(os.path.join(imgDir, "pickledDetections.pkl"), "rb")

    detectedFiberList = pickle.load(pkl_file)
    pkl_file.close()

    # firstFrameNum = frameNumFromName(detectedFiberList[0]["imageFrames"][0])
    # lastFrameNum = frameNumFromName(detectedFiberList[-1]["imageFrames"][-1])
    nDetectedFibers = len(detectedFiberList)
    print("found", nDetectedFibers, "fibers")

    # # create visualization
    # imgArray = numpy.zeros((nDetectedFibers, lastFrameNum-firstFrameNum+1))
    # fig = plt.figure(figsize=(20,5))
    # for fiberNum, detectedFiber in enumerate(detectedFiberList):
    #     for counts, filename in itertools.izip(detectedFiber["counts"], detectedFiber["imageFrames"]):
    #         frameNumber = frameNumFromName(filename)
    #         imgInd = frameNumber - firstFrameNum
    #         imgArray[fiberNum, imgInd] = counts / float(numpy.max(detectedFiber["counts"]))
    #         # plt.plot(fiberNum, imgInd, ".k")

    # plt.imshow(imgArray, interpolation="none")
    # nfn = os.path.join(imgDir, "fiberMap.png")
    # plt.show()

    # fig = plt.figure(figsize=(20,5))
    # for fiber in imgArray:
    #     plt.plot(range(lastFrameNum-firstFrameNum+1), fiber, 'k')
    # plt.show()

    SlitheadSolver(detectedFiberList)

    # maxDetect, imgFrames = extractMax(detectedFiberList)
    # plt.hist(numpy.diff(imgFrames), bins=200)
    # plt.show()

    # findGaps(detectedFiberList)
    # import pdb; pdb.set_trace()
    #fig.savefig(nfn)
    #plt.close(fig)    # close the figure
