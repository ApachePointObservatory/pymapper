from __future__ import division, absolute_import

import itertools
import copy
import os

import numpy

import scipy.ndimage
import logging

from astropy.io import fits

from pymapper.camera import FScanCamera, pickleDetectionList, unpickleDetectionList, MINCOUNTS, MINSEP, DetectedFiber, unpickleCentroids
from pymapper.fiberAssign import SlitheadSolver, FocalSurfaceSolver

# from sdss.utilities.yanny import yanny

class FScanFrame(object):
    def __init__(self, frame, motor, motorpos, tstamp, row, col, flux, flatImg=None):
        """
        Class representing a single fscan frame
        @param[in] frame: int, frame number
        @param[in] flatImg: 2D array. A flat image from the CCD.
        @param[in] motor: int, motor ID (1-3)
        @param[in] tstamp: float, timestamp for frame
        @param[in] row: int, row of pixel
        @param[in] col: int, col of pixel
        @param[in] flux: int, flux of pixel
        """
        self.frame = frame
        if flatImg is not None:
            self.flatImg = flatImg
        else:
            self.flatImg = numpy.ones((960,960))*5
        self.motor = motor
        self.motorpos = motorpos / 34.2857 + 40.44 # scale factor from APO motor steps to mm
        self.tstamp = tstamp
        self.rows = [row]
        self.cols = [col]
        self.fluxes = [flux]

    def addPixel(self, row, col, flux):
        """Add a single pixel and it's flux

        @param[in] row: int, row of pixel
        @param[in] col: int, col of pixel
        @param[in] flux: int, flux of pixel
        """
        self.rows.append(row)
        self.cols.append(col)
        self.fluxes.append(flux)

    def getImg(self):
        """Return a 2D image of this frame
        Add fluxes from specfic pixels to the flat image

        from IDL:
       ; New camera installed 02-Jul-2009; scalex=-0.754, scaley=-0.752 mm/pix
       camparam_new = $
        {scalex: [-0.78,-0.72,0.005], $ ; range of X scales
         xpixfac: -1.0              , $ ; factor for xpix (invert X axis)
         nmotor: 65000L             , $ ; max motor steps + buffer for x-correlation
         msigma: 18.                  } ; Gaussian width of spot in motor steps

        """
        imgData = copy.deepcopy(self.flatImg) + numpy.random.standard_normal((960,960))
        for row, col, flux in itertools.izip(self.rows, self.cols, self.fluxes):
            # imgData[col*-1, row] += flux
            imgData[row, col] += flux
        # if self.frame in [3939]:
        #     hdu = fits.PrimaryHDU(imgData)
        #     hdulist = fits.HDUList([hdu])
        #     filename = "fits%i.fits"%self.frame
        #     if os.path.exists(filename):
        #         os.remove(filename)
        #     hdulist.writeto(filename)
        return imgData

class ScanMovie(object):
    def __init__(self, imgDir, fscanFile, flatImgFile=None):
        """Class for representing an evil mapper fscan as a pymapper movie

        @param[in]  fscanFile: path to fscan file
        @param[in]  flatImgFile: path to a flatImgFile
        """
        self.imgDir = imgDir
        self.fscanFile = fscanFile
        if flatImgFile is not None:
            self.flatImg = scipy.ndimage.imread(flatImgFile)
        else:
            self.flatImg = numpy.ones((960,960))*5

        self.frames = self.framesFromFile(fscanFile)
        # self._flatFile = None


    # @property
    # def flatFile(self):
    #     if self._flatFile is None:
    #         flatFile = numpy.zeros((960,960))
    #         for ii, bias in enumerate(self.frames[200:400]):
    #             flatFile = flatFile + bias.getImg()
    #         flatFile / float(ii+1)
    #         self._flatFile = flatFile
    #     return self._flatFile

    def parseLine(self, line):
        """Parse a single line beginning with SCANPIX
        """
        spix, motor, frame, motorpos, tstamp, row, col, flux = line.split()
        return (int(motor), int(frame), int(motorpos),
            float(tstamp), int(row), int(col), int(flux))

    def framesFromFile(self, fscanFile, stack=1):
        """Parse the fscan file return an array of frame objects
        """
        frames = []
        with open(fscanFile) as f:
            filelines = f.readlines()
        currentFrame = None
        currentFrameNum = 0
        for line in filelines:
            if not line.startswith("SCANPIX"):
                continue
            motor, frame, motorpos, tstamp, row, col, flux = self.parseLine(line)
            # print("frame: ", frame)
            if motor == 1:
                continue # next line until we hit motor 2 which is apogee
            if motor == 3:
                break # don't continue reading lines
            # only want motor == 2 (apogee motor)
            # if frame > 700:
            # #     # print("coninuing", frame, motor)
            #     break
            if currentFrame is None:
                # first frame
                currentFrameNum = frame
                currentFrame = FScanFrame(frame, motor, motorpos, tstamp, row, col, flux, self.flatImg)
            elif frame == currentFrame.frame*stack + stack:
                # new frame
                frames.append(currentFrame)
                currentFrameNum += 1
                print("new frame at %i"%frame)
                currentFrame = FScanFrame(currentFrameNum, motor, motorpos, tstamp, row, col, flux, self.flatImg)
            else:
                # append to current frame
                currentFrame.addPixel(row, col, flux)
        # add final frame
        # might miss one frame?
        frames.append(currentFrame)
        return frames

def sortDetections(fScanCamera, plot=False, minCounts=MINCOUNTS, minSep=MINSEP):
    """Reorganize detection list into groups
    of detections (1 group per fiber)
    """
    detectedFibers = []
    for brightestCentroid in fScanCamera.centroidList:
        isNewDetection = None
        crashMe = False
        if brightestCentroid["counts"] is not None and brightestCentroid["counts"] > MINCOUNTS:
            isNewDetection = True
            # search through every previous detection
            for prevDetection in detectedFibers:
                if prevDetection.belongs2me(brightestCentroid, minSep):
                    if isNewDetection == False:
                        crashMe = True
                        logging.warn("bad bad, crash me! in %s"%brightestCentroid["imageFile"])
                    # print("previous detection!!!", brightestCentroid.xyCtr)
                    isNewDetection = False
                    # print("previous detection", os.path.split(imageFile)[-1], prevDetection.imageFiles)
                    prevDetection.add2me(brightestCentroid)
                    # break # assign to the first that works???
                # was this is a new detection?
            if isNewDetection:
                # print('new detection:', os.path.split(imageFile)[-1], brightestCentroid.counts, brightestCentroid.xyCtr)
                detectedFibers.append(DetectedFiber(brightestCentroid))

        if plot:
            frameNumber = brightestCentroid["frame"]
            print("plotting frame %i"%frameNumber)
            color = "r" if isNewDetection else "b"
            fig = plt.figure(figsize=(10,10));plt.imshow(fScanCamera.frames[frameNumber].getImg(), vmin=0, vmax=10)#plt.show(block=False)
            plt.scatter(0, 0, s=80, facecolors='none', edgecolors='b')
            if brightestCentroid["xyCtr"] is not None:
               x,y = brightestCentroid["xyCtr"]
               plt.scatter(x, y, s=80, facecolors='none', edgecolors=color)
            zfilled = "%i"%frameNumber
            zfilled = zfilled.zfill(5)
            dd = os.path.split(fScanCamera.imageDir)[0]
            nfn = os.path.join(dd, "pyguide%s.png"%zfilled)
            fig.savefig(nfn); plt.close(fig)    # close the figure
        # if crashMe:
        #     raise RuntimeError("Non-unique detection!!!!")
    return detectedFibers

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    pyMapperDir = os.getenv("PYMAPPER_DIR")
    testDir = os.path.join(pyMapperDir, "tests/data/")
    fscanFile = os.path.join(testDir, "8637/fiberScan-8637-57534-01.par")
    flatFile = os.path.join(testDir, "flat.bmp")
    plPlugMapFile = os.path.join(testDir, "8637", "plPlugMapM-8637-57534-01.par")
    sm = ScanMovie(testDir, fscanFile, flatImgFile=None)
    fsc = FScanCamera(testDir, sm)
    # fsc.reprocessImages()
    fsc.frames = sm.frames
    # detectedFiberList = sortDetections(fsc, plot=False)
        # pickle and save the detection list
    # pickleDetectionList(detectedFiberList, testDir)
    plt.imshow(sm.frames[0].getImg())
    # plt.show()
    rawfluxes = []
    rawMotorPos = []
    countsList = []
    centroidList = unpickleCentroids(testDir)
    detectionList = unpickleDetectionList(testDir)
    fig = plt.figure()
    for centroid, frame in itertools.izip(centroidList, sm.frames):
        rawflux = numpy.sum(frame.fluxes)
        counts = centroid["counts"]
        if counts > 0:
            counts = rawflux
        rawfluxes.append(rawflux)
        countsList.append(counts)
        rawMotorPos.append(frame.motorpos)
    imgNums = []
    detectionCounts = []
    detectionMotorPos = []
    nFrameList = []
    for detection in detectionList:
        nFrames = len(detection.imageFiles)
        nFrameList.append(nFrames)
        frameInts = [int(x.split(".")[0]) for x in detection.imageFiles]
        middleFrame = nFrames // 2
        detectionCounts.append(detection.counts[middleFrame]+15000)
        imgNums.append(frameInts[middleFrame])
        detectionMotorPos.append(detection.motorPos)

    # plt.plot(range(len(countsList)), countsList, '+r') #countsList ~ frames not mm so doesn't scale right anymore
    # plt.plot(range(len(rawfluxes)), rawfluxes)
    plt.plot(rawMotorPos, rawfluxes)
    plt.plot(detectionMotorPos, detectionCounts, 'or')
    print(len(detectionList), "detections")

    # next determine slithead solution:
    APOfiberpos = os.path.join(os.getenv("PYMAPPER_DIR"), "etc", "fiberslitposAPO.par")
    shs = SlitheadSolver(detectionList, APOfiberpos)
    shs.getOffsetAndScale()
    shs.matchDetections()
    fss = FocalSurfaceSolver(detectionList, plPlugMapFile)
    # centroidList = unpickleCentroids(testDir)
    # fluxes = []
    # fig = plt.figure()
    # for centroid in centroidList:
    #     fluxes.append(centroid["counts"])
    # plt.plot(range(len(fluxes)), fluxes)

    plt.show()

    import pdb; pdb.set_trace()