""" camera control
"""
from __future__ import division, absolute_import

import os
import glob
import subprocess
import pickle
import time
import logging
import traceback
from multiprocessing import Pool
import multiprocessing as mp

import shutil

import scipy.ndimage
import numpy
from astropy.io import fits


import PyGuide

from twisted.internet import reactor

from . import plt

import cProfile, pstats, StringIO

# how to deal with hot pixels?
"""
from IDL:
   ; Search for hot pixels, and subtract the median value (plus some
   ; multiple of the dispersion of that value) from those pixels
   ; A pixel is hot if it's detected in >5% of the timesteps
   nframetot = n_elements(uniq(sdat.frame))
   xypix = sdat.col * fscanRows + sdat.row
   isort = sort(xypix)
   iuniq = uniq(xypix[isort])
   nuniq = n_elements(iuniq)
   i0 = [0,iuniq[0:nuniq-1]+1]
   num = iuniq - i0 + 1
   ihot = where(num GT 0.05*nframetot, nhot)
   for i=0L, nhot-1L do begin
      thisxy = xypix[isort[i0[ihot[i]]]]
      thisindx = isort[i0[ihot[i]]:iuniq[ihot[i]]]
      djs_iterstat, sdat[thisindx].flux, sigrej=3, $
       median=thismed, sigma=thissig
      hotval = thismed + 3*thissig + 5
      sdat[thisindx].flux = (sdat[thisindx].flux - hotval) > 0
      splog, 'Hot pixel at X,Y=', thisxy/fscanRows, thisxy MOD fscanRows, $
       ' val=', hotval
   endfor

   ; Invert camera image if need be
   sdat.col *= camparam.xpixfac

"""

# from pymapper.imgProcess import multiprocessImage
# how to ensure that the image program is killed?

# EXE = "testFakeWrite.py"
# BASENAME = "pyTestFile"
# IMGEXTENSION = "foo"

EXE = "AsynchronousGrabWrite"
IMGBASENAME = "img"
IMGEXTENSION = "fits"
MINCOUNTS = 0
ROUGH_THRESH = 8 # for filtering out pixels going into "totalCounts" in each frame
MINSEP = 3.5 # min between fibers separation in pixels

# CCDInfo = PyGuide.CCDInfo(bias=50, readNoise=10, ccdGain=1)
CCDInfo = PyGuide.CCDInfo(bias=1, readNoise=0, ccdGain=1)
# pyguide mask
PyGuideMask = numpy.zeros((960,960))
midPix = 960/2.
for x in range(960):
    for y in range(960):
        pixRad = numpy.sqrt((x-midPix)**2+(y-midPix)**2)
        if pixRad > midPix:
            PyGuideMask[x,y]=1

# self.imageDir = os.path.join(os.environ["PYMAPPER_DIR"], "tests")

# globals requred due to multiprocessing breaking
# when calling class methods rather than functions
# accepting a single argument (processImage)
# TZERO = None # time stamp of first image
MOTORSPEED = None
MOTORSTART = None
MOTOREND = None

def getScanParams(paramfile):
    """Parse the paramfile to determine the scan params
    """
    outDict = {
        "speed": None,
        "start": None,
        "end": None,
        "plateID": None,
    }
    with open(paramfile, "r") as f:
        paramLine = f.readlines()

    for line in paramLine:
        param, val = line.split()
        if param == "plateID":
            outDict["plateID"] = int(val)
        elif param == "startPos":
            outDict["start"] = float(val)
        elif param == "endPos":
            outDict["end"] = float(val)
        elif param == "speed":
            outDict["speed"] = float(val)
    if None in outDict.values():
        raise RuntimeError("Could not extract plateID, start, end, and/or speed from paramfile")
    return outDict

def _basePickle(basename, pyobj, scanDir):
    # save detections to picked file
    filename = os.path.join(scanDir, "%s.pkl"%basename)
    if os.path.exists(filename):
        # rename existing pickle file with a timestamp
        # when it was last modified
        fTime = time.ctime(os.path.getmtime(filename))
        movedfilename = "%s-%s.pkl"%(basename, fTime)
        shutil.move(filename, os.path.join(scanDir, movedfilename))
        print("moving %s to %s"%(filename, movedfilename))
    output = open(filename, "wb")
    pickle.dump(pyobj, output)
    print("dumping to %s"%(filename))
    output.close()

def _baseUnpickle(basename, scanDir):
   pkl = open(os.path.join(scanDir, "%s.pkl"%basename), "rb")
   pyObj = pickle.load(pkl)
   pkl.close()
   return pyObj

def pickleCentroids(centroidList, scanDir):
    # save detections to picked file
    return _basePickle("centroidList", centroidList, scanDir)

def unpickleCentroids(scanDir):
    return _baseUnpickle("centroidList", scanDir)

def pickleDetectionList(detectionList, scanDir):
    return _basePickle("detectionList", detectionList, scanDir)

def unpickleDetectionList(scanDir):
    return _baseUnpickle("detectionList", scanDir)

# def frameNumFromName(imgName, imgBase=IMGBASENAME, imgExt=IMGEXTENSION):
#     imgName = os.path.split(imgName)[-1]
#     return int(imgName.split(imgBase)[-1].split(".%s"%imgExt)[0])

# def convToFits(imageFileDirectory, flatImg, frameStartNum, frameEndNum=None, imgBaseName=IMGBASENAME, imgExtension=IMGEXTENSION):
#     saveImage()

# def saveImage(filename, array2d):
#     if os.path.exists(filename):
#         print("removing previous", filename)
#         os.remove(filename)
#     # median filter the image
#     # array2d = scipy.ndimage.median_filter(array2d, size=1)
#     if filename.endswith(".fits"):
#         _saveFits(filename, array2d)
#     else:
#         _saveJpeg(filename, array2d)

# def _saveFits(filename, array2d):
#     hdu = fits.PrimaryHDU(array2d)
#     hdu.writeto(filename)

# def _saveJpeg(filename, array2d):
#     scipy.misc.imsave(filename, array2d)

class Camera(object):
    def __init__(self, imageDir, motorStart, motorEnd, motorSpeed):
        global MOTORSTART
        global MOTOREND
        global MOTORSPEED
        MOTORSTART = motorStart
        MOTOREND = motorEnd
        MOTORSPEED = motorSpeed
        self.acquiring = False
        self.process = None
        self.imageDir = imageDir
        # assert os.path.exists(imageDir), "%s doesn't exist, create it"%imageDir
        # whine if this directory already has images
        # assert len(self.getAllImgFiles())==0, "%s is not empty!"%imageDir
        self.acquisionCB = None
        self.procImgCall = None
        self.centroidList = []

    def doneProcessingCallback(self, callFunc):
        # to be called when all image processing is done,
        # the centroids should have been saved as a picklefile
        self.procImgCall = callFunc

    def getAllImgFiles(self):
        # glob doesn't order correctly in
        return sorted(glob.glob(os.path.join(self.imageDir, "*."+IMGEXTENSION)))

    def getNthFile(self, fileNum):
        filenumStr = ("%i"%fileNum).zfill(6)
        filename = "%s%s.%s"%(IMGBASENAME, filenumStr, IMGEXTENSION)
        return os.path.join(self.imageDir, filename)

    def getUnprocessedFileList(self):
        # order matters!
        return [self.getNthFile(fileNum) for fileNum in range(self.nFilesProcessed+1, self.nFilesWritten+1)]

    @property
    def nFilesProcessed(self):
        return len(self.centroidList)

    @property
    def nFilesWritten(self):
        return len(self.getAllImgFiles())

    def beginAcquisition(self, callFunc=None):
        """call callFunc when acquision has started (first image seen in directory)
        """
        print("beginAcquisition")
        if callFunc is not None:
            print("setting acquision cb: %s"%str(callFunc))
            self.acquisitionCB = callFunc
        self.acquiring = True
        # this process initializes and starts the camera
        self.process = subprocess.Popen(EXE, cwd=self.imageDir)
        self.waitForFirstImage()

    def stopAcquisition(self):
        print("Stopping Camera Acquision")
        self.process.kill()
        self.acquiring = False

    def multiprocessDone(self):
        print("All frames processed!")
        print("pickling centroid list")
        totalTime = time.time()-self.processingStart
        print("total time: %.2f"%totalTime)
        print("processed %.2f frames per second"% (len(self.centroidList)/totalTime))
        pickleCentroids(self.centroidList, self.imageDir)
        if self.procImgCall is not None:
            self.procImgCall()

    def waitForFirstImage(self):
        # loop here until the first image is seen.  Once it is
        # fire the acquisition callback and begin processing
        # global TZERO
        if os.path.exists(self.getNthFile(1)):
            print("acquisition started")
            # set the zeroth timestamp for determining motor position for each image
            # TZERO = os.path.getmtime(self.getNthFile(1))
            if self.acquisitionCB is not None:
                print("firing acquisition callback")
                reactor.callLater(0., self.acquisitionCB)
            # the first image is here, we're free to start
            # processing them
            reactor.callLater(0., self.multiprocessImageLoop)
            self.processingStart = time.time()
        else:
            # first image not seen yet try again
            reactor.callLater(0., self.waitForFirstImage)


    def reprocessImages(self):
        # global TZERO
        # self.tZero = os.path.getmtime(self.getNthFile(1))
        # TZERO = self.tZero
        self.processingStart = time.time()
        allImgFiles = self.getAllImgFiles()
        self.centroidList = multiprocessImage(allImgFiles, callFunc=None, block=True)
        self.multiprocessDone()

    def multiprocessImageLoop(self, centroidList=None):
        # called recursively until all images are (multi!) processed
        print("multiprocessImageLoop")
        if centroidList:
            print("adding %i centroids"%len(centroidList))
            self.centroidList.extend(centroidList)
        unprocessedFileList = self.getUnprocessedFileList()
        # don't process more than 50 images at a time
        unprocessedFileList = unprocessedFileList[:50]
        if unprocessedFileList:
            print("processing images %s to %s"%tuple([os.path.split(_img)[-1] for _img in [unprocessedFileList[0], unprocessedFileList[-1]]]))
            multiprocessImage(unprocessedFileList, self.multiprocessImageLoop, block=False)
        else:
            # no files to process.
            print("no files to process")
            if self.acquiring:
                print("still acquiring")
                # camera is still acquiring, so continue calling myself
                reactor.callLater(0., self.multiprocessImageLoop)
            else:
                # camera is done, no remaining files to process
                self.multiprocessDone()


def multiprocessImage(imageFileList, callFunc=None, block=False):
    # may want to try map_async
    p = Pool(5)
    if block:
        output = p.map(processImage, imageFileList)
        if callFunc is not None:
            callFunc(output)
        return output
    else:
        return p.map_async(processImage, imageFileList, callback=callFunc)

def processImage(imageFile):
    """! Process a single image

    @param[in] imageFile. String

    """
    # global TZERO
    global MOTORSPEED
    global MOTORSTART
    global MOTOREND
    motorDirection = numpy.sign(MOTOREND-MOTORSTART)
    # frame = int(imageFile.split(IMGBASENAME)[-1].split(".")[0])
    # timestamp = os.path.getmtime(imageFile) - TZERO
    frame = None
    timestamp = None
    counts = None
    xyCtr = None
    rad = None
    totalCounts = None

    # put some filtering here
    # try:

    fitsImg = fits.open(imageFile)
    timestamp = fitsImg[0].header["TSTAMP"]
    frame = fitsImg[0].header["FNUM"]
    imgData = fitsImg[0].data
    fitsImg.close()

    if frame % 100 == 0:
        print("processing frame number %i"%frame)
    flatImg = imgData.flatten()
    thresh = flatImg[numpy.nonzero(flatImg>ROUGH_THRESH)]
    # print("median %.4f thresh %.4f  %i pixels over thresh"%(medianValue, sigma, len(thresh)))
    totalCounts = numpy.sum(thresh)

    pyGuideCentroids = PyGuide.findStars(imgData, PyGuideMask, None, CCDInfo)[0]
    # did we get any centroids?
    if pyGuideCentroids:
        counts = pyGuideCentroids[0].counts
        xyCtr = pyGuideCentroids[0].xyCtr
        rad = pyGuideCentroids[0].rad
    del imgData # paranoia

    # except Exception:
    #     print("some issue with pyguide on img (skipping): ", imageFile)
    #     traceback.print_exc()

    return dict((
                    ("imageFile", imageFile),
                    ("counts", counts),
                    ("xyCtr", xyCtr),
                    ("rad", rad),
                    ("motorPos", MOTORSTART + motorDirection*MOTORSPEED*timestamp),
                    ("frame", frame),
                    ("totalCounts", totalCounts),
                ))

def processImage_(args):
    out = [None]
    cProfile.runctx('out[0] = processImage(args)', globals(), locals(),
                    'profile-%s.out' % mp.current_process().name)
    return out[0]

class FScanCamera(Camera):
    """For testing with existing idlmapper fcan files
    """
    def __init__(self, imageDir, scanMovie):
        """@param[in] scanMovie a ScanMovie Obj
        """
        self.imageDir = imageDir
        self.scanMovie = scanMovie


    def reprocessImages(self):
        self.centroidList = self.multiprocessImage(self.scanMovie.frames, callFunc=None, block=True)
        pickleCentroids(self.centroidList, self.imageDir)

    def multiprocessImage(self, fScanFrames, callFunc, block=False):
        # may want to try map_async
        p = Pool(2)
        if block:
            output = p.map(self.processImage, fScanFrames)
            return output
        else:
            return p.map_async(self.processImage, fScanFrames, callback=callFunc)

    def processImage(self, fScanFrame):
        """! Process a single image

        @param[in] fScanFrame: an FScanFrame obj

        """
        # print("processing img: ", os.path.split(imageFile)[-1])
        if fScanFrame.frame % 100 == 0:
            print("frame %i   %.2f done"%(fScanFrame.frame, float(fScanFrame.frame+1)/float(len(self.scanMovie.frames))*100))
        imgData = fScanFrame.getImg()
        # if fScanFrame.frame > 0:
        #     imgData = imgData - self.scanMovie.frames[fScanFrame.frame-1].getImg()#/ self.scanMovie.flatFile
        imageFile = "%i.fscanframe"%fScanFrame.frame
        motorPos = fScanFrame.motorpos
        counts = None
        xyCtr = None
        rad = None
        totalCounts = numpy.sum(imgData)
        # k = numpy.array([[.1,.1,.1],[.1,1,.1],[.1,.25,.1]])
        # imgData = scipy.ndimage.convolve(imgData, k)
        # imgData = scipy.ndimage.filters.median_filter(imgData, 3)
        try:
            pyGuideCentroids = PyGuide.findStars(imgData, None, None, CCDInfo)[0]
            # did we get any centroids?
            if pyGuideCentroids:
                counts = pyGuideCentroids[0].counts
                xyCtr = pyGuideCentroids[0].xyCtr
                rad = pyGuideCentroids[0].rad
        except Exception:
            print("some issue with pyguide on img (skipping): ", imageFile)
            traceback.print_exc()
        return dict((
                        ("imageFile", imageFile),
                        ("counts", counts),
                        ("xyCtr", xyCtr),
                        ("rad", rad),
                        ("motorPos", motorPos),
                        ("frame", fScanFrame.frame),
                        ("totalCounts", totalCounts)
                    ))

class DetectedFiber(object):
    def __init__(self, centroidDict):
        self.centroidList = [centroidDict]
        # self.centroids = [pyGuideCentroid]
        # self.imageFiles = [imageFileName]
        self.slitHeadInd = None
        self.plPlugObjInd = None
        self.focalPos = [None, None]

    def setSlitIndex(self, index):
        self.slitHeadInd = index

    def getSlitIndex(self):
        return self.slitHeadInd

    def setPlPlugObjInd(self, index):
        self.plPlugObjInd = index

    def getPlPlugObjInd(self):
        return self.plPlugObjInd

    def setFocalPos(self, focalPos):
        self.focalPos = focalPos

    @property
    def imageFiles(self):
        return [centroid["imageFile"] for centroid in self.centroidList]

    @property
    def counts(self):
        return [centroid["counts"] for centroid in self.centroidList]

    @property
    def xyCtr(self):
        # return center based on weighted counts
        return numpy.average([cent["xyCtr"] for cent in self.centroidList], axis=0, weights=self.counts)

    @property
    def xyCtrs(self):
        return [centroid["xyCtr"] for centroid in self.centroidList]

    @property
    def motorPos(self):
        # return center based on weighted counts
        return numpy.average([cent["motorPos"] for cent in self.centroidList], axis=0, weights=self.counts)

    @property
    def motorPositions(self):
        return [centroid["motorPos"] for centroid in self.centroidList]

    @property
    def frame(self):
        # return center based on weighted counts
        return numpy.average([cent["frame"] for cent in self.centroidList], axis=0, weights=self.counts)

    @property
    def frames(self):
        return [centroid["frame"] for centroid in self.centroidList]

    @property
    def rad(self):
        return numpy.average([cent["rad"] for cent in self.centroidList], axis=0, weights=self.counts)

    def belongs2me(self, centroidDict, minSep=MINSEP):
        # if center moves by more than 0.25 pixels
        # doesn't belong
        dist = numpy.linalg.norm(numpy.subtract(centroidDict["xyCtr"], self.xyCtr))
        # if dist < 3:
        #     print("belongs to", dist, self.imageFiles)
        return dist < minSep
        # print("dist!", dist, self.imageFiles)
        # return dist<(self.rad/2.)

    def add2me(self, centroidDict):
        self.centroidList.append(centroidDict)

    def detectedIn(self, imageFileName):
        return imageFileName in self.imageFiles


def sortDetections(brightestCentroidList, plot=False, minCounts=MINCOUNTS, minSep=MINSEP):
    """Reorganize detection list into groups
    of detections (1 group per fiber)
    """
    detectedFibers = []
    prevFrame = -1
    # get the motor positions
    timestamps = [cent["motorPos"] for cent in brightestCentroidList]
    tsdiff = numpy.diff(timestamps)
    print("motorpos/frame: %.4f +/- %.4f"%(numpy.mean(tsdiff), numpy.std(tsdiff)))
    for brightestCentroid in brightestCentroidList:
        if brightestCentroid["frame"] < prevFrame:
            raise RuntimeError("centroids not sorted in ascending frame order!")
        prevFrame = brightestCentroid["frame"]
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
            imageFile = brightestCentroid["imageFile"]
            color = "r" if isNewDetection else "b"
            img = fits.open(imageFile)
            imgData = img[0].data
            frameNumber = img[0].header["FNUM"]
            img.close()
            fig = plt.figure(figsize=(10,10));plt.imshow(imgData, vmin=0, vmax=10)#plt.show(block=False)
            plt.scatter(0, 0, s=80, facecolors='none', edgecolors='b')
            if brightestCentroid["xyCtr"] is not None:
               x,y = brightestCentroid["xyCtr"]
               plt.scatter(x, y, s=80, facecolors='none', edgecolors=color)
            zfilled = "%i"%frameNumber
            zfilled = zfilled.zfill(5)
            dd = os.path.split(imageFile)[0]
            nfn = os.path.join(dd, "pyguide%s.png"%zfilled)
            fig.savefig(nfn); plt.close(fig)    # close the figure
        # if crashMe:
        #     raise RuntimeError("Non-unique detection!!!!")
    # throw out any false detections (must show up in at least two consecutive frames)
    filteredDetectedFibers = []
    for detection in detectedFibers:
        if len(detection.frames) < 3:
            print("Removing false detecion in frames:", detection.frames)
            continue
        filteredDetectedFibers.append(detection)
    return filteredDetectedFibers

if __name__ == "__main__":
    camera = Camera()
    camera.start()
    reactor.callLater(3*60., camera.stop)
    reactor.run()