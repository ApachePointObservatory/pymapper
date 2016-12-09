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

from . import plt

# from .measureSlitPos import getMeasuredFiberPositions

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
        # self.hack()
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

    def correlate(self):
        maxCor = None
        maxShift = None
        maxScale = None
        rescaledMotorPositions = None
        # for scale in numpy.arange(1-0.002, 1+0.002, 0.0005):
        for scale in numpy.arange(0.9, 1.1, 0.0005):
            # print("scale: ", scale)
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
        self.scale = maxScale
        print("scale: %.8f offset: %.4f"%(self.scale, self.offset))
        # self.rescaledMotorPositions = rescaledMotorPositions - self.offset
        self.rescaledMotorPositions = numpy.asarray([cent["motorPos"] for cent in self.centroidList])*self.scale - self.offset
        self.scaledDetections = numpy.asarray([det.motorPos for det in self.detectedFiberList])*self.scale - self.offset
        print("got %i detections"%len(self.scaledDetections))

    def hack(self):
        self.nomMotorPos = numpy.asarray(self.nomMotorPos) * 1.05 #- 4

    def plotSolution(self, scanDir=None):
        # scale rawFluxes
        # rawFluxes = (self.rawFluxes - numpy.median(self.rawFluxes))
        # scale by 3/4 the max value
        scaleValue = sorted(self.rawFluxes)[int(len(self.rawFluxes)*3/4)]
        rawFluxes = 0.5*self.rawFluxes / scaleValue
        modelMotorPos = numpy.arange(self.rescaledMotorPositions[0], self.rescaledMotorPositions[-1], -1*FIBERDIAMETER/30)
        unscaledDetections = numpy.asarray([det.motorPos for det in self.detectedFiberList])
        fiberNums = [det.getSlitIndex() for det in self.detectedFiberList]
        modeledFlux = [self.fluxFromMotorPos(motorPos) for motorPos in modelMotorPos]
        fig = plt.figure(figsize=(200,30))
        # import pdb; pdb.set_trace()
        slitModel, = plt.plot(modelMotorPos, modeledFlux, '-k', alpha=0.5, linewidth=3)
        scaledRaw, = plt.plot(self.rescaledMotorPositions, rawFluxes, '.-g', alpha=0.5, linewidth=3)
        scaledDetections, = plt.plot(self.rescaledMotorPositions, self.normalizedFlux, '.-b', alpha=0.5, linewidth=3)
        for unscaled, scaled in itertools.izip(unscaledDetections, self.scaledDetections):
            # plt.plot(unscaled, 0.6, 'xr')
            # plt.plot(scaled, 0.51, 'or')
            plt.plot([unscaled, scaled], [.6, .51], 'r-', alpha=0.5, linewidth=2)
        unscaledCenters, = plt.plot(unscaledDetections, 0.6*numpy.ones(len(unscaledDetections)), 'xr')
        scaledCenters, = plt.plot(self.scaledDetections, 0.51*numpy.ones(len(self.scaledDetections)), 'or')
        for x, fiberNum in itertools.izip(self.scaledDetections, fiberNums):
            plt.text(x, 0.52, str(fiberNum+1), horizontalalignment='center', fontsize=8)
        if scanDir:
            nfn = os.path.join(scanDir, "slitheadSolution.png")
        else:
            nfn = "slitheadSolution.png"
        plt.legend(
            [slitModel, scaledRaw, scaledDetections, unscaledCenters, scaledCenters],
            ["Gaussian Slit Model", "Fit, Thresholded, Normed Total Counts In Frame", "Fit, Normed Detection Counts", "Unfit Detection Centers", "Fit Detection Centers"],
            fontsize=30,
            # loc = "upper center",
            )
        plt.xlabel("Motor Position (mm)")
        plt.ylabel("Normalized Counts")
        meanX = numpy.mean(self.rescaledMotorPositions)
        maxY = numpy.max(rawFluxes)
        plt.text(meanX, maxY - 0.1, "RMS after scale/offeset fitting: %.6f"%self.rms, horizontalalignment="center", fontsize=30)
        missingFiberStr = ",".join([str(fiber) for fiber in self.missingFibers])
        plt.text(meanX, maxY - 0.15, "Missing Fiber Numbers: %s"%missingFiberStr, horizontalalignment="center", fontsize=30)
        plt.text(meanX, maxY - 0.05, "Fit -- Scale: %.4f Offset: %.4f (mm)"%(self.scale, self.offset), horizontalalignment="center", fontsize=30)
        fig.savefig(nfn); plt.close(fig)
        # plt.show()


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
        normalizedFlux = numpy.zeros(self.centroidList[-1]["frame"]-1)
        rawFlux = numpy.array([cent["totalCounts"] for cent in self.centroidList])
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
        # measMatchPos = [detectedFiber.motorPos/self.scale - self.offset for detectedFiber in self.detectedFiberList]
        inds = []
        errs = []
        for measMotorPos in self.scaledDetections:
            # for each measured position find the closest match
            # to a nominal fiber (after scale and offset have been applied)
            # it is assued that the measured positions should be
            # very close to only one of the nominal positions
            allErrs = numpy.abs(measMotorPos-self.nomMotorPos)
            minInd = numpy.argmin(allErrs)
            err = allErrs[minInd]
            errs.append(err)
            assert minInd not in inds # can have on fiber match twice
            if len(inds)>0:
                assert minInd > inds[-1]
            inds.append(minInd)
        self.matchInds = inds
        self.missingFibers = sorted([fiber+1 for fiber in set(range(300))-set(self.matchInds)])
        # next calculate the RMS from all errors
        self.rms = numpy.sqrt(numpy.sum(numpy.asarray(errs)**2) / len(errs)) # rms err in mm
        # set each detected fibers pos on the slit head
        for matchInd, detectedFiber in itertools.izip(self.matchInds, self.detectedFiberList):
            detectedFiber.setSlitIndex(matchInd)




class PlPlugMap(object):
    def __init__(self, plPlugMapFile, scanDir, cartID, fscanID, fscanMJD):
        self.scanDir = scanDir
        self.cartID = cartID
        self.fscanID = fscanID
        self.fscanMJD = fscanMJD
        self.plPlugMap = yanny(filename=plPlugMapFile, np=True)
        self.objectInds = numpy.argwhere(self.plPlugMap["PLUGMAPOBJ"]["holeType"]=="OBJECT")
        self.xPos = self.plPlugMap["PLUGMAPOBJ"]["xFocal"][self.objectInds].flatten()
        self.yPos = self.plPlugMap["PLUGMAPOBJ"]["yFocal"][self.objectInds].flatten()
        self.radPos = numpy.sqrt(self.xPos**2+self.yPos**2)
        self.plateID = int(self.plPlugMap["plateId"])

    def append(self, updateDict):
        self.plPlugMap.append(updateDict)

    def enterMappedData(self, detectedFiberList):
        # for fibers not found enter fiberID = -1, spectrographID = -1, and throughput = 0
        # first set spectrograph id and fiber id to -1 for all OBJECTS
        # self.plPlugMap["PLUGMAPOBJ"]["spectrographId"][self.objectInds] = 6
        # self.plPlugMap["PLUGMAPOBJ"]["fiberId"][self.objectInds] = -1
        # throughput should already be 0....do we need to check?\
        for detectedFiber in detectedFiberList:
            objInd = detectedFiber.getPlPlugObjInd()
            plInd = self.objectInds[objInd]
            self.plPlugMap["PLUGMAPOBJ"]["fiberId"][plInd] = detectedFiber.getSlitIndex()
            self.plPlugMap["PLUGMAPOBJ"]["throughput"][plInd] = detectedFiber.counts
            # paranoia!
            assert self.plPlugMap["PLUGMAPOBJ"][plInd]["holeType"].flatten()[0] == "OBJECT"

    # def _enterMappedData(self, mappedFiberInds, mappedFiberThroughputs, mappedPlPlugObjInds):
    #     # for fibers not found enter fiberID = -1, spectrographID = -1, and throughput = 0
    #     # first set spectrograph id and fiber id to -1 for all OBJECTS
    #     # self.plPlugMap["PLUGMAPOBJ"]["spectrographId"][self.objectInds] = 6
    #     # self.plPlugMap["PLUGMAPOBJ"]["fiberId"][self.objectInds] = -1
    #     # throughput should already be 0....do we need to check?\
    #     for plInd, fiberInd, fiberThroughput in itertools.izip(mappedPlPlugObjInds, mappedFiberInds, mappedFiberThroughputs):
    #         # fiberInd is fiber number + 1 (zero indexed)
    #         # plInd is index based on holeType==OBJECT only (we threw away all other types)
    #         # first find out what plInd index corresponds to in scope of all holes in file
    #         holeInd = self.objectInds[plInd]
    #         self.plPlugMap["PLUGMAPOBJ"]["spectrographId"][holeInd] = 0
    #         # because of zero index add 1, because 1-16 are reserved for guide fibers begin count at 17
    #         self.plPlugMap["PLUGMAPOBJ"]["fiberId"][holeInd] = fiberInd + 1 + 16
    #         self.plPlugMap["PLUGMAPOBJ"]["throughput"][holeInd] = fiberThroughput

    def writeMe(self):
        # replace plPlugMapP-XXX with plPlugMapM-XXX
        # previousDirectory, filename = os.path.split(self.plPlugMap.filename)
        # determine any existing scans
        fscanIDstr = ("%i"%self.fscanID).zfill(2)
        fileName = "plPlugMapM-%i-%i-%s.par"%(self.plateID, self.fscanMJD, fscanIDstr)
        updateDict = {
            "fscanMJD": self.fscanMJD,
            "fscanId": self.fscanID,
            "cartridgeId": self.cartID,
            "instruments": "APOGEE_SOUTH"
            }
        self.plPlugMap.append(updateDict)
        self.plPlugMap.write(os.path.join(self.scanDir, fileName))

    def plotMissing(self):
        tabHeight = 0.5 * MMPERINCH
        tabWidth = 0.8 * MMPERINCH
        fig = plt.figure(figsize=(50,50))
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
            if fiberID < 0:
                # fiber not found
                marker = "x" # red x
                color = "red"
            else:
                # fiber found
                marker = "o" #circle
                color = "black"
            plt.plot(xPos, yPos, marker=marker, color=color, fillstyle="none", mew=2)
            # plt.text(xPos*-1, yPos, str(objInd[0]+1), horizontalalignment="center", fontsize=10)
        # plt.show(block=True)
        figname = os.path.join(self.scanDir, "unplugged.png")
        fig.savefig(figname); plt.close(fig)



class FocalSurfaceSolver(object):
    def __init__(self, detectedFiberList, plPlugMapFile, scanDir, cartID, fscanID, fscanMJD):
        self.scanDir = scanDir
        self.detectedFiberList = detectedFiberList
        if not os.path.exists(plPlugMapFile):
            raise RuntimeError("coudl not locate plPlugMapFile: %s"%plPlugMapFile)
        self.plPlugMap = PlPlugMap(plPlugMapFile, scanDir, cartID, fscanID, fscanMJD)
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
        self.plPlugMap.enterMappedData(self.detectedFiberList)
        self.plPlugMap.writeMe()
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
        # measXPos = xPos - numpy.mean(xPos)
        measXPos = -1*(xPos - numpy.mean(xPos))
        # measYPos = -1*(yPos - numpy.mean(yPos))
        measYPos = yPos - numpy.mean(yPos)
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

    def nextPlotFilename(self):
        baseName = "transrotscale"
        ii = 0
        while True:
            zfillNum = str(ii).zfill(2)
            filename = os.path.join(self.scanDir, baseName+"%s.png"%zfillNum)
            if os.path.exists(filename):
                ii += 1
                continue
            else:
                break
        return filename


    def plot(self, measx, measy, block=False):
        fig = plt.figure(figsize=(30,30))
        plt.plot(self.plPlugMap.xPos, self.plPlugMap.yPos, 'or', alpha=0.5)
        plt.plot(measx, measy, 'ob', alpha=0.5)
        plt.xlim([-400, 400])
        plt.ylim([-400, 400])
        fig.savefig(self.nextPlotFilename())
        plt.show(block=block)
        plt.close(fig)

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
                print("got a multimatch")
                multiMatchInds.append(plInd)
                badIndex = plPlugMapInds.index(plInd)
                plPlugMapInds.pop(badIndex)
                measPosInds.pop(badIndex)
            else:
                plPlugMapInds.append(plInd)
                measPosInds.append(measInd)

        # did we get the expected amount of matches?
        print("got ", len(plPlugMapInds), "matches", len(multiMatchInds), "multimatches")
        if currentCall>1 and len(plPlugMapInds) == len(xArray) or currentCall == maxCalls:
            # must iterate at least once to fully solve using russells model!
            # every match found, we're done!
            # self.plPlugMapInds = plPlugMapInds
            # self.measPosInds = measPosInds
            self.measXPos = xArray
            self.measYpos = yArray
            errArr = []
            for measPosInd, plInd in itertools.izip(measPosInds, plPlugMapInds):
                detectedFiber = self.detectedFiberList[measPosInd]
                x = xArray[measPosInd]
                y = yArray[measPosInd]
                detectedFiber.setPlPlugObjInd(plInd)
                detectedFiber.setFocalPos([x,y])
                xPl = self.plPlugMap.xPos[plInd]
                yPl = self.plPlugMap.yPos[plInd]
                xyErr = (x-xPl)**2+(y-yPl)**2
                # print("xyErr %.4f"%xyErr)
                errArr.append(xyErr)
            plt.hist(errArr, 100)
            print("XY RMS: %.4f"%(numpy.sqrt(numpy.sum(errArr)/len(errArr))))
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
        # if previous solution is same as new solution ==  no change, try to eliminate false positives?
        trans, rot, scale = transRotScaleSolution
        if trans[0] == trans[1] == rot == 0 and scale == 1:
            # find n(multimatch) number of points furthest from any neighbor
            # and remove them
            print("no change in transrotscale, removing multimatches")
            distList = [0]*len(multiMatchInds)
            maxInds = [-1]*len(multiMatchInds)
            for measInd, (xPos, yPos) in enumerate(itertools.izip(xArray, yArray)):
                # copied from above, break out?
                dist = numpy.min(numpy.sqrt((self.plPlugMap.xPos - xPos)**2+(self.plPlugMap.yPos - yPos)**2))
                if dist > numpy.max(distList):
                    minInd = numpy.argmin(distList)
                    distList[minInd] = dist
                    maxInds[minInd] = measInd
            for measInd, dist in itertools.izip(maxInds, distList):
                print("removing detection %i with nearest neighbor dist %.4f"%(measInd, dist))
            raise RuntimeError("Implement me!?")

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
