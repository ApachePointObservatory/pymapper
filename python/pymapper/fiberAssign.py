"""Routines for identifying, assigning fibers. Determining which (if any are missing)


todo:
in recursive fitter, identify when solution is not changing and exit
what to do with multiple matches?
discover and reject when too bright?rawImage-8787-57476-shortexp bad, _slow_dark ok
set timeouts!!!1
"""
from __future__ import division, absolute_import

import os
import itertools
import glob
import logging

import numpy
from scipy.optimize import fmin
from scipy.optimize import minimize_scalar
import scipy.optimize


from sdss.utilities.yanny import yanny
from fitPlugPlateMeas.fitData import TransRotScaleModel, ModelFit

# from .measureSlitPos import getMeasuredFiberPositions

# spacing between blocks (measured empirically), units are adjacent fiber spacing
# turns out to be better to just use the same spacing for every block...
# BlockSpaceMultiplier = [1.928, 1.868, 1.828, 1.877, 1.910, 1.805, 1.980, 1.917, 1.887]

DEBUG = False
# fwhm are used to build gaussians in a minimizing "energy" function
# for fitting scale, tranlation, rotation to fiber positions
FWHMCOARSE = 5 # mm
FWHMFINE = 0.1 # mm
MMPERINCH = 25.4 # mm per inch
MICRONSPERMM = 1000.0
# plates have 32 inch diameter
PLATERADIUS = 32 * MMPERINCH / 2.
# MATCHTHRESHb = 3 #mm (require a match to be within this threshold to be considered robust)
FIBERDIAMETER = 120 / MICRONSPERMM
SLIT_MATCH_THRESH = FIBERDIAMETER / 2. #mm wiggle room for slit head matching 1/2 of a fiber diameter
INTER_SLIT_FIBERCENTER_DIST = 350 / MICRONSPERMM
INTER_BLOCK_FIBERCETNER_DIST = 0.5 #mm

class SlitheadSolver(object):
    def __init__(self, detectedFiberList, centroidList, fiberslitposFile=None):
        """List of detections
        """
        if fiberslitposFile is None:
            fiberslitposFile = os.path.join(os.getenv("PYMAPPER_DIR"), "etc", "fiberslitpos.dat")
        self.detectedFiberList = detectedFiberList
        self.centroidList = centroidList
        self.nomFiberNum, self.nomMotorPos = self.parseFiberPosFile(fiberslitposFile)
        self.hackNomPos()
        # det motor pos is motor position at each frame
        # normailzed flux is counts in each detection normalized by
        # total counts in the detection (so the sum of counts in any detection is 1)
        # normalized flux only contians counts that belong to a detection
        # raw fluxes is just the sum of counts in each frame.
        self.detMotorPos, self.normalizedFlux, self.rawFluxes = self.generateDetectionTrace()
        # self.modeledMotorPos = numpy.arange(self.nomMotorPos[0]-4, self.nomMotorPos[-1]+4, FIBERDIAMETER/4)
        # self.modeledFluxes = [self.fluxFromMotorPos(motorPos) for motorPos in self.modeledMotorPos]

        # rescaled motor positions are set in correlate
        # they are used for the actual matching of detected fibers to the model
        # of fibers on the slit. rescaled motor positions are
        # self.detMotorPos with the scale and offset applied
        self.rescaledMotorPositions = None
        self.correlate()

    def hackNomPos(self):
        self.nomMotorPos = numpy.asarray(self.nomMotorPos) * 1.05 - 1 + 10

    def correlate(self):
        maxCor = None
        maxShift = None
        maxScale = None
        rescaledMotorPositions = None
        for scale in numpy.arange(1-0.002, 1+0.002, 0.0005):
            modeledMotorPos = self.detMotorPos * scale
            modeledFluxes = [self.fluxFromMotorPos(motorPos) for motorPos in modeledMotorPos]
            cout = numpy.correlate(self.normalizedFlux, modeledFluxes, mode="full")
            maxInd = numpy.argmax(cout)
            maxVal = cout[maxInd]
            shift = maxInd - int(len(cout)/2) # shift is in units of frames
            if maxCor is None or maxVal > maxCor:
                maxCor = maxVal
                maxShift = shift
                maxScale = scale
                rescaledMotorPositions = modeledMotorPos
        # print("shift of %i for scale %.8f"%(maxShift, maxScale))
        self.offset = maxShift * numpy.mean(numpy.diff(self.detMotorPos*maxScale)) # convert frames to motor pos
        print("offset: %.4f"%self.offset)
        self.scale = maxScale
        self.rescaledMotorPositions = rescaledMotorPositions + self.offset


    def plotSolution(self):


    def parseFiberPosFile(self, fiberslitposFile):
        fiberNums = []
        motorPositions = []
        with open(fiberslitposFile, "r") as f:
            filelines = f.readlines()
        for line in filelines:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fiberNum, motorPos = line.split()
            fiberNums.append(int(fiberNum))
            motorPositions.append(float(motorPos))
        return fiberNums, motorPositions

    def generateDetectionTrace(self):
        normalizedFlux = numpy.zeros(self.centroidList[-1]["frame"])
        rawFlux = numpy.arrray([cent["totalCounts"] for cent in self.centroidList])
        detMotorPos = numpy.array([cent["motorPos"] for cent in self.centroidList])
        for detection in self.detectedFiberList:
            # normalize all counts to sum to 1
            # for each detection
            detSum = numpy.sum(detection.counts)
            for frame, counts in itertools.izip(detection.frames, detection.counts):
                normalizedFlux[frame-1] = counts/detSum
        return detMotorPos, normalizedFlux, rawFlux


    def fluxFromMotorPos(self, motorPos):
        """Return expected power from a motorPosition based on
        the measured positions on the slit
        model as a gaussain at each measured fiber position on the slit
        """
        minDist = numpy.min(numpy.abs(numpy.asarray(self.nomMotorPos) - motorPos))
        # print("minDist", minDist)
        # thats the minimum distance to the nearest gaussian (or expected position)
        return 0.5*numpy.exp(-1*(minDist)**2/float(FIBERDIAMETER**2))


    def matchDetections(self):
        # scale the measured motor positions of by the scaling
        # nomMatchPos = numpy.asarray(self.nomMotorPos)*self.scale + self.offset
        measMatchPos = [detectedFiber.motorPos*self.scale + self.offset for detectedFiber in self.detectedFiberList]
        inds = []
        errs = []
        for measMotorPos in measMatchPos:
            # for each measured position find the closest match
            # to a nominal fiber (after scale and offset have been applied)
            # it is assued that the measured positions should be
            # very close to only one of the nominal positions
            allErrs = numpy.abs(measMotorPos-self.rescaledMotorPositions)
            minInd = numpy.argmin(allErrs)
            err = allErrs[minInd]
            errs.append(err)
            if minInd in inds:
                import pdb; pdb.set_trace()
            assert minInd not in inds # can have on fiber match twice
            inds.append(minInd)
        self.matchInds = inds
        self.missingFibers = sorted([fiber+1 for fiber in set(range(300))-set(self.matchInds)])
        # next calculate the RMS from all errors
        self.rms = numpy.sqrt(numpy.sum(numpy.asarray(errs)**2) / len(errs)) # rms err in mm
        # set each detected fibers pos on the slit head
        for matchInd, detectedFiber in itertools.izip(self.matchInds, self.detectedFiberList):
            detectedFiber.setSlitIndex(matchInd)




class PlPlugMap(object):
    def __init__(self, plPlugMapFile, scanDir):
        self.scanDir = scanDir
        self.plPlugMap = yanny(filename=plPlugMapFile, np=True)
        self.objectInds = numpy.argwhere(self.plPlugMap["PLUGMAPOBJ"]["holeType"]=="OBJECT")
        self.xPos = self.plPlugMap["PLUGMAPOBJ"]["xFocal"][self.objectInds].flatten()
        self.yPos = self.plPlugMap["PLUGMAPOBJ"]["yFocal"][self.objectInds].flatten()
        self.radPos = numpy.sqrt(self.xPos**2+self.yPos**2)
        self.plateID = int(self.plPlugMap["plateId"])

    def enterMappedData(self, detectedFiberList):
        # for fibers not found enter fiberID = -1, spectrographID = -1, and throughput = 0
        # first set spectrograph id and fiber id to -1 for all OBJECTS
        # self.plPlugMap["PLUGMAPOBJ"]["spectrographId"][self.objectInds] = 6
        # self.plPlugMap["PLUGMAPOBJ"]["fiberId"][self.objectInds] = -1
        # throughput should already be 0....do we need to check?\
        for detectedFiber in detectedFiberList:
            objInd = detectedFiber.getPlPlugObjInd()
            plInd = self.objectInds[plInd]
            self.plPlugMap["PLUGMAPOBJ"]["fiberId"][plInd] = detectedFiber.getSlitIndex()
            self.plPlugMap["PLUGMAPOBJ"]["throughput"][plInd] = detectedFiber.counts

    def _enterMappedData(self, mappedFiberInds, mappedFiberThroughputs, mappedPlPlugObjInds):
        # for fibers not found enter fiberID = -1, spectrographID = -1, and throughput = 0
        # first set spectrograph id and fiber id to -1 for all OBJECTS
        # self.plPlugMap["PLUGMAPOBJ"]["spectrographId"][self.objectInds] = 6
        # self.plPlugMap["PLUGMAPOBJ"]["fiberId"][self.objectInds] = -1
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
        # print("globStr", globStr)
        nExisting = glob.glob(globStr)
        # number this scan accordingly
        scanNum = len(nExisting) + 1
        scanNumStr = ("%i"%scanNum).zfill(2)
        # construct the new file name, will not overwrite any existing (scanNumStr is incremented)
        filename = "plPlugMapM-%i-%i-%s-test.par"%(self.plateID, mjd, scanNumStr)
        self.plPlugMap.write(os.path.join(writeDir, filename))

    def plotMissing(self):
        tabHeight = 0.5 * MMPERINCH
        tabWidth = 0.8 * MMPERINCH
        fig = plt.figure(figsize=(10,10))
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
        # plt.show(block=True)
        figname = os.path.join(self.scanDir, "unplugged.png")
        fig.savefig(figname); plt.close(fig)



class FocalSurfaceSolver(object):
    def __init__(self, detectedFiberList, plPlugMapFile, scanDir):
        self.scanDir = scanDir
        self.detectedFiberList = detectedFiberList
        if not os.path.exists(plPlugMapFile):
            raise RuntimeError("coudl not locate plPlugMapFile: %s"%plPlugMapFile)
        self.plPlugMap = PlPlugMap(plPlugMapFile, scanDir)
        #self.nHistBins = 100 # number of bins for histogram
        #self.focalRadHist, self.binEdges = numpy.histogram(self.plPlugMap.radPos, bins=self.nHistBins)
        # self.plotRadialHist(self.plPlugMap.radPos, bins=self.binEdges)
        # compute radial histogram of x y focal positions
        # self.plotRadialHist(measR, bins=self.binEdges)
        self.measXPos, self.measYPos = self.initialMeasTransforms()
        self.fitTransRotScale()
        self.matchMeasToPlPlugMap(self.measXPos, self.measYPos) # sets attriubte plPlugMapInds
        # throughputList = self.getThroughputList()
        # self.plPlugMap.enterMappedData(self.measPosInds, throughputList, self.plPlugMapInds)
        self.plPlugMap.enterMappedData(self.detectedFiberList, self.plPlugMapInds)
        self.plPlugMap.writeMe(scanDir, 55555)
        self.plPlugMap.plotMissing()

    def getThroughputList(self):
        # report max counts as throughput, what should the true definition be?
        return [numpy.max(self.detectedFiberList[ind].counts) for ind in self.measPosInds]

    def initialMeasTransforms(self):
        # just put the measured data in the ballpark of the
        # focal positions
        # for initial rough scale and x,y tranlation fitting to
        # measured data
        # return x,y coordinate list for detections
        # "centroid" detection by weighted counts..
        detectionCenters = numpy.asarray([detection.xyCtr for detection in self.detectedFiberList])
        # for detectedFiber in self.detectedFiberList:
        #     xyCenters = [numpy.asarray(xyCenter) for xyCenter in detectedFiber.xyCtrs]
        #     counts = detectedFiber.counts
        #     detectionCenters.append(numpy.average(xyCenters, axis=0, weights=counts))
        # detectionCenters = numpy.asarray(detectionCenters)
        xPos, yPos = detectionCenters[:,0], detectionCenters[:,1]
        # pyguide convention +y goes down the image.  this is bad, flip it here
        # for matching to the plPlugMap focal values
        # also determine the rough plate center, averaging x and y
        measXPos = xPos - numpy.mean(xPos)
        measYPos = -1*(yPos - numpy.mean(yPos))
        # measYPos = measYPos - numpy.mean(measYPos)
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

    # def getMeasuredCenters(self):
    #     # return x,y coordinate list for detections
    #     # "centroid" detection by weighted counts..
    #     detectionCenters = numpy.asarray([detection.xyCtr for detection in self.detectedFiberList])
    #     # for detectedFiber in self.detectedFiberList:
    #     #     xyCenters = [numpy.asarray(xyCenter) for xyCenter in detectedFiber.xyCtrs]
    #     #     counts = detectedFiber.counts
    #     #     detectionCenters.append(numpy.average(xyCenters, axis=0, weights=counts))
    #     # detectionCenters = numpy.asarray(detectionCenters)
    #     xPos, yPos = detectionCenters[:,0], detectionCenters[:,1]
    #     return xPos, yPos

    def matchMeasToPlPlugMap(self, xArray, yArray, currentCall=0, maxCalls=10, previousSolution=None):
        """Match measured positions to
        those in the plPlugMap file.  Exactly (in a least squares sense) solve for
        translation rotation and scale using only robust matches.  Iterate until we have all of em.

        recursively!!!, this routine tweaks trans, rot, and scale each time
        """
        if currentCall == maxCalls:
            # raise RuntimeError("Max recursion reached!!!!")
            print("max recurion reached")
        #     return
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
        if len(plPlugMapInds) == len(xArray) or currentCall == maxCalls:
            # every match found, we're done!
            self.plPlugMapInds = plPlugMapInds
            # self.measPosInds = measPosInds
            self.measXPos = xArray
            self.measYpos = yArray
            for detectedFiber,x,y,plInd in itertools.izip(self.detectedFiberList, xArray, yArray, plPlugMapInds):
                detectedFiber.setPlPlugObjInd(plInd)
                detectedFiber.setFocalPos([x,y])
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
        self.plot(xArray, yArray)
        self.matchMeasToPlPlugMap(xArray, yArray, currentCall=currentCall+1, previousSolution=transRotScaleSolution)



        # plt.figure()
        # plt.hist(rads, bins=100)
        # plt.show(block=False)

        # plt.figure()
        # plt.hist(err, bins=100)
        # plt.show(block=True)

if __name__ == "__main__":
    import pickle
    imgDir = "/Volumes/Boof/scan/57476/rawImage-8787-57476-shortexp_slow_dark"# "/Users/csayres/Desktop/mapperPyguide/run6noFlat"

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

    # SlitheadSolver(detectedFiberList)
    fss = FocalSurfaceSolver(detectedFiberList, os.path.join(imgDir, "plPlugMapP-8787.par"))
    # maxDetect, imgFrames = extractMax(detectedFiberList)
    # plt.hist(numpy.diff(imgFrames), bins=200)
    # plt.show()

    # findGaps(detectedFiberList)
    # import pdb; pdb.set_trace()
    #fig.savefig(nfn)
    #plt.close(fig)    # close the figure
