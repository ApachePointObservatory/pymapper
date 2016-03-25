"""Routines for identifying, assigning fibers. Determining which (if any are missing)
"""
from __future__ import division, absolute_import

import os
import itertools
import glob

import numpy
from scipy.optimize import fmin

# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sdss.utilities.yanny import yanny
from fitPlugPlateMeas.fitData import TransRotScaleModel, ModelFit

# spacing between blocks (measured empirically), units are adjacent fiber spacing
# turns out to be better to just use the same spacing for every block...
# BlockSpaceMultiplier = [1.928, 1.868, 1.828, 1.877, 1.910, 1.805, 1.980, 1.917, 1.887]

DEBUG = False
# fwhm are used to build gaussians in a minimizing "energy" function
# for fitting scale, tranlation, rotation to fiber positions
FWHMCOARSE = 5 # mm
FWHMFINE = 0.1 # mm
MMPERINCH = 25.4 # mm per inch
# plates have 32 inch diameter
PLATERADIUS = 32 * MMPERINCH / 2.
# MATCHTHRESH = 3 #mm (require a match to be within this threshold to be considered robust)

class ModeledTrace(object):
    def __init__(self, fiberSpacing, blockSpacing):
        """Model for a slithead trace, units are frames for both fiber and blockSpacing
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

class SlitheadSolver(object):
    def __init__(self, detectedFiberList):
        """List of detections
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
        for detectedFiber in detectedFiberList:
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
                        # i suppose this could happen if pixels are saturated
                        normalizedCounts = 0.999
                    else:
                        foundMax = True
                assert measuredTrace[traceInd] == 0, "already counts found in this frame?!?!"
                measuredTrace[traceInd] = normalizedCounts
        return measuredTrace

class PlPlugMap(object):
    def __init__(self, plPlugMapFile):
        self.plPlugMap = yanny(filename=plPlugMapFile, np=True)
        self.objectInds = numpy.argwhere(self.plPlugMap["PLUGMAPOBJ"]["holeType"]=="OBJECT")
        self.xPos = self.plPlugMap["PLUGMAPOBJ"]["xFocal"][self.objectInds].flatten()
        self.yPos = self.plPlugMap["PLUGMAPOBJ"]["yFocal"][self.objectInds].flatten()
        self.radPos = numpy.sqrt(self.xPos**2+self.yPos**2)
        self.plateID = int(self.plPlugMap["plateId"])

    def enterMappedData(self, mappedFiberInds, mappedFiberThroughputs, mappedPlPlugObjInds):
        # for fibers not found enter fiberID = -1, spectrographID = -1, and throughput = 0
        # first set spectrograph id and fiber id to -1 for all OBJECTS
        self.plPlugMap["PLUGMAPOBJ"]["spectrographId"][self.objectInds] = -1
        self.plPlugMap["PLUGMAPOBJ"]["fiberId"][self.objectInds] = -1
        # throughput should already be 0....do we need to check?\
        for plInd, fiberInd, fiberThroughput in itertools.izip(mappedPlPlugObjInds, mappedFiberInds, mappedFiberThroughputs):
            # fiberInd is fiber number + 1 (zero indexed)
            # plInd is index based on holeType==OBJECT only (we threw away all other types)
            # first find out what plInd index corresponds to in scope of all holes in file
            holeInd = self.objectInds[plInd]
            self.plPlugMap["PLUGMAPOBJ"]["spectrographId"][holeInd] = 0
            # because of zero index add 1, because 1-16 are reserved for guide fibers begin count at 17
            self.plPlugMap["PLUGMAPOBJ"]["fiberId"][holeInd] = fiberInd + 1 + 16
            self.plPlugMap["PLUGMAPOBJ"]["throughput"][holeInd] = fiberThroughput

    def writeMe(self, writeDir, mjd):
        # replace plPlugMapP-XXX with plPlugMapM-XXX
        # previousDirectory, filename = os.path.split(self.plPlugMap.filename)
        # determine any existing scans
        globStr = os.path.join(writeDir, "plPlugMapM-%i-%i-*.par"%(self.plateID, mjd))
        print("globStr", globStr)
        nExisting = glob.glob(globStr)
        # number this scan accordingly
        scanNum = len(nExisting) + 1
        scanNumStr = ("%i"%scanNum).zfill(2)
        # construct the new file name, will not overwrite any existing (scanNumStr is incremented)
        filename = "plPlugMapM-%i-%i-%s.par"%(self.plateID, mjd, scanNumStr)
        self.plPlugMap.write(os.path.join(writeDir, filename))

    def plotMissing(self):
        tabHeight = 0.5 * MMPERINCH
        tabWidth = 0.8 * MMPERINCH
        plt.figure(figsize=(10,10))
        plt.box(on=True)
        radLimit = PLATERADIUS + tabHeight + 10 #mm
        limits = (-1*radLimit, radLimit)
        plt.xlim(limits)
        plt.ylim(limits)
        # plot the plate circle
        thetas = numpy.linspace(0, 2*numpy.pi, 1000)
        plateX = PLATERADIUS * numpy.cos(thetas)
        plateY = PLATERADIUS * numpy.sin(thetas)
        # determine angle from tabWitdh
        totalCircumference = 2 * numpy.pi * PLATERADIUS
        tabFraction = tabWidth / totalCircumference
        tabAngularFraction = tabFraction * numpy.pi * 0.95 # by eye tweaking
        # replace the radii around the tab with something longer
        tabArgs = numpy.nonzero(numpy.logical_and(thetas < 1.5 * numpy.pi + tabAngularFraction, thetas > 1.5 * numpy.pi - tabAngularFraction))
        plateX[tabArgs] = (PLATERADIUS+tabHeight)*numpy.cos(thetas[tabArgs])
        plateY[tabArgs] = (PLATERADIUS+tabHeight)*numpy.sin(thetas[tabArgs])
        # diameter is 32 inches
        plt.plot(plateX, plateY, "-k")
        # plot found fibers:
        for objInd in self.objectInds:
            xPos = self.plPlugMap["PLUGMAPOBJ"]["xFocal"][objInd]
            yPos = self.plPlugMap["PLUGMAPOBJ"]["yFocal"][objInd]
            fiberID = self.plPlugMap["PLUGMAPOBJ"]["fiberId"][objInd]
            if fiberID == -1:
                # fiber not found
                marker = "x" # red x
                color = "red"
            else:
                # fiber found
                marker = "o" #circle
                color = "black"
            plt.plot(xPos, yPos, marker=marker, color=color, fillstyle="none", mew=2)
        plt.show(block=True)



class FocalSurfaceSolver(object):
    def __init__(self, detectedFiberList, plPlugMapFile):
        self.detectedFiberList = detectedFiberList
        self.plPlugMap = PlPlugMap(plPlugMapFile)
        #self.nHistBins = 100 # number of bins for histogram
        #self.focalRadHist, self.binEdges = numpy.histogram(self.plPlugMap.radPos, bins=self.nHistBins)
        # self.plotRadialHist(self.plPlugMap.radPos, bins=self.binEdges)
        # compute radial histogram of x y focal positions
        # self.plotRadialHist(measR, bins=self.binEdges)
        self.measXPos, self.measYPos = self.initialMeasTransforms()
        self.fitTransRotScale()
        self.matchMeasToPlPlugMap(self.measXPos, self.measYPos) # sets attriubte plPlugMapInds
        throughputList = self.getThroughputList()
        self.plPlugMap.enterMappedData(self.measPosInds, throughputList, self.plPlugMapInds)
        self.plPlugMap.writeMe("/Users/csayres/Desktop/", 55555)
        self.plPlugMap.plotMissing()

    def getThroughputList(self):
        # report max counts as throughput, what should the true definition be?
        return [numpy.max(self.detectedFiberList[ind]["counts"]) for ind in self.measPosInds]

    def initialMeasTransforms(self):
        # just put the measured data in the ballpark of the
        # focal positions
        # for initial rough scale and x,y tranlation fitting to
        # measured data
        measXPos, measYPos = self.getMeasuredCenters()
        # pyguide convention +y goes down the image.  this is bad, flip it here
        # for matching to the plPlugMap focal values
        # also determine the rough plate center, averaging x and y
        measXPos = measXPos - numpy.mean(measXPos)
        measYPos = -1*(measYPos - numpy.mean(measYPos))
        # in polar coords..
        measR = numpy.sqrt(measXPos**2+measYPos**2)
        measTheta = numpy.arctan2(measYPos, measXPos)
        # determine a rough scaling based on x,y value range in plPlugMap
        # file
        maxPlRadPos = numpy.max(self.plPlugMap.radPos)
        maxMeasRadPos = numpy.max(measR)
        roughScale = maxPlRadPos / maxMeasRadPos
        # apply the scale to the measured (polar) positions
        measR = measR * roughScale
        measXPos = measR * numpy.cos(measTheta)
        measYPos = measR * numpy.sin(measTheta)
        return measXPos, measYPos

    def plot(self, measx, measy, block=False):
        plt.figure()
        plt.plot(self.plPlugMap.xPos, self.plPlugMap.yPos, 'or')
        plt.plot(measx, measy, 'ok')
        plt.show(block=block)

    def fitTransRotScale(self):
        """Use a minimizer to get translation, rotation and scale close enough
        to determine (hopefully many) robust matches.
        """
        self.plot(self.measXPos, self.measYPos)
        # step 1 rough fit tranlation
        transx, transy = fmin(self.minimizeTranslation, [0,0], args=(FWHMCOARSE,))
        self.plot(self.measXPos-transx, self.measYPos-transy)

        # step 2 rough fig scale and rot
        rot, scale = fmin(self.minimizeRotScale, [0,1], args=(transx, transy, FWHMCOARSE))
        x, y = self.applyTransRotScale(self.measXPos, self.measYPos, transx, transy, rot, scale)
        self.plot(x, y)
        print(transx, transy, rot, scale)

        # step 3, re fit trans rot scale together, with tighter gaussians around the target points
        transx, transy, rot, scale = fmin(self.minimizeTransRotScale, [transx, transy, rot, scale], args=(FWHMFINE,))
        print(transx, transy, rot, scale)
        self.measXPos, self.measYPos = self.applyTransRotScale(self.measXPos, self.measYPos, transx, transy, rot, scale)
        self.plot(self.measXPos, self.measYPos)

    def minimizeTranslation(self, transxy, fwhm):
        transx, transy = transxy
        x, y = self.applyTransRotScale(self.measXPos, self.measYPos, transx, transy)
        return self.computeEnergy(x, y, fwhm)

    def minimizeRotScale(self, rotScale, transx, transy, fwhm):
        rot, scale = rotScale
        x, y = self.applyTransRotScale(self.measXPos, self.measYPos, transx, transy, rot, scale)
        return self.computeEnergy(x, y, fwhm)

    def minimizeTransRotScale(self, transRotScale, fwhm):
        transx, transy, rot, scale = transRotScale
        x, y = self.applyTransRotScale(self.measXPos, self.measYPos, transx, transy, rot, scale)
        return self.computeEnergy(x, y, fwhm)

    def computeEnergy(self, xPos, yPos, fwhm = 5):
        """place gaussian over each point in plplugmap file
        fwhm is mm (units of plplugmap)
        """
        # compute all pairwise differences
        # could be vectorized if needed...
        # after all this function is evaluated a lot by the minimizer.
        totalEnergy = 0
        for x, y in itertools.izip(xPos, yPos):
            dist = numpy.sqrt((self.plPlugMap.xPos-x)**2+(self.plPlugMap.yPos-y)**2)
            # energy is a gaussian, lower distance == higher energy
            energy = numpy.sum(numpy.exp(-1 * dist**2/(2*fwhm**2)))
            totalEnergy += energy
        # return negative energy because we are minimizing this
        # print totalEnergy
        return -1 * totalEnergy

    def applyTransRotScale(self, x, y, transx, transy, rot=0, scale=1):
        x = x - transx
        y = y - transy
        if rot != 0 or scale != 1:
            r = numpy.sqrt(x**2+y**2)
            theta = numpy.arctan2(y, x)
            theta = theta - rot
            r = r * scale
            x = r * numpy.cos(theta+2*numpy.pi)
            y = r * numpy.sin(theta+2*numpy.pi)
        return x, y

    def getMeasuredCenters(self):
        # return x,y coordinate list for detections
        # "centroid" detection by weighted counts..
        detectionCenters = []
        for detectedFiber in self.detectedFiberList:
            xyCenters = [numpy.asarray(xyCenter) for xyCenter in detectedFiber["xyCtrs"]]
            counts = detectedFiber["counts"]
            detectionCenters.append(numpy.average(xyCenters, axis=0, weights=counts))
        detectionCenters = numpy.asarray(detectionCenters)
        xPos, yPos = detectionCenters[:,0], detectionCenters[:,1]
        return xPos, yPos

    def matchMeasToPlPlugMap(self, xArray, yArray, currentCall=0, maxCalls=10, previousSolution=None):
        """Match measured positions to
        those in the plPlugMap file.  Exactly (in a least squares sense) solve for
        translation rotation and scale using only robust matches.  Iterate until we have all of em.

        recursively!!!, this routine tweaks trans, rot, and scale each time
        """
        if currentCall == maxCalls:
            raise RuntimeError("Max recursion reached!!!!")
        # brute force, compare distance to
        # all other points, could
        # do some nearest neighbor search
        # or could vecorize
        multiMatchInds = []
        plPlugMapInds = []
        measPosInds = []
        for measInd, (xPos, yPos) in enumerate(itertools.izip(xArray, yArray)):
            dist = numpy.sqrt((self.plPlugMap.xPos - xPos)**2+(self.plPlugMap.yPos - yPos)**2)
            plInd = numpy.argmin(dist)
            if plInd in plPlugMapInds:
                # this index was already matched to another
                # put it in the multimatch index
                multiMatchInds.append(plInd)
                badIndex = plPlugMapInds.index(plInd)
                plPlugMapInds.pop(badIndex)
                measPosInds.pop(badIndex)
            else:
                plPlugMapInds.append(plInd)
                measPosInds.append(measInd)

        # did we get the expected amount of matches?
        print("got ", len(plPlugMapInds), "matches", len(multiMatchInds), "multimatches")
        if len(plPlugMapInds) == len(xArray):
            # every match found, we're done!
            self.plPlugMapInds = plPlugMapInds
            self.measPosInds = measPosInds
            self.measXPos = xArray
            self.measYpos = yArray
            return

        # determine exact transrotscale solution now that we have (at least) some robust matches.
        measPos = numpy.asarray([ [xArray[ii], yArray[ii]] for ii in measPosInds ])
        nomPos = numpy.asarray([ [self.plPlugMap.xPos[ii], self.plPlugMap.yPos[ii]] for ii in plPlugMapInds ])
        fitTransRotScale = ModelFit(
            model = TransRotScaleModel(),
            measPos = measPos,
            nomPos = nomPos,
            doRaise=True,
        )
        transRotScaleSolution = fitTransRotScale.model.getTransRotScale()
        # is this the same solution as the previous? if so complain now (we aren't getting anywhere)
        # @todo, implement later!!!!
        print("iter", currentCall, "fit model", transRotScaleSolution)
        # apply this new (tweaked) model to every point and call this routine again
        # (until we get all matches found)
        xyArray = numpy.asarray([xArray, yArray]).T
        newPositions = fitTransRotScale.model.apply(xyArray, doInverse=True)
        xArray = newPositions[:,0]
        yArray = newPositions[:,1]
        self.matchMeasToPlPlugMap(xArray, yArray, currentCall=currentCall+1, previousSolution=transRotScaleSolution)



        # plt.figure()
        # plt.hist(rads, bins=100)
        # plt.show(block=False)

        # plt.figure()
        # plt.hist(err, bins=100)
        # plt.show(block=True)

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
    fss = FocalSurfaceSolver(detectedFiberList, os.path.join(imgDir, "plPlugMapP-8787.par"))
    # maxDetect, imgFrames = extractMax(detectedFiberList)
    # plt.hist(numpy.diff(imgFrames), bins=200)
    # plt.show()

    # findGaps(detectedFiberList)
    # import pdb; pdb.set_trace()
    #fig.savefig(nfn)
    #plt.close(fig)    # close the figure
