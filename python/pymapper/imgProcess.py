"""Routines for processing images from the lco mapper
"""
from __future__ import division, absolute_import
import glob
import os
import time
from multiprocessing import Pool
import traceback

import numpy
import scipy.ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.io import fits

import PyGuide

THRESH = 50
MINCOUNTS = 0
MINSEP = 3.5 # min between fibers separation in pixels

CCDInfo = PyGuide.CCDInfo(bias=50, readNoise=10, ccdGain=1)

def processImage(imageFile):
    """! Process a single image

    @param[in] imageFile. String

    """
    # logging.info("processing img: ", os.path.split(imageFile)[-1])
    imgData = scipy.ndimage.imread(imageFile)
    counts = None
    xyCtr = None
    rad = None
    try:
        pyGuideCentroids = PyGuide.findStars(imgData, None, None, CCDInfo)[0]
        # did we get any centroids?
        if pyGuideCentroids:
            counts = pyGuideCentroids[0].counts
            xyCtr = pyGuideCentroids[0].xyCtr
            rad = pyGuideCentroids[0].rad
    except Exception as e:
        logging.info("some issue with pyguide on img (skipping): ", imageFile)
        traceback.logging.info_exc()
    return dict((
                    ("imageFile", imageFile),
                    ("counts", counts),
                    ("xyCtr", xyCtr),
                    ("rad", rad)
                ))

def applyThreshold(array2d, thresh):
    # which pixels are greater than sigmaDetect sigma from the mean
    pixInd = numpy.argwhere(array2d < thresh)
    array2d[pixInd] = 0
    return array2d

def to2d(array1d, imshape):
    return numpy.reshape(array1d, imshape)

def saveImage(filename, array2d):
    if os.path.exists(filename):
        logging.info("removing previous", filename)
        os.remove(filename)
    # median filter the image
    # array2d = scipy.ndimage.median_filter(array2d, size=1)
    if filename.endswith(".fits"):
        _saveFits(filename, array2d)
    else:
        _saveJpeg(filename, array2d)

def _saveFits(filename, array2d):
    hdu = fits.PrimaryHDU(array2d)
    hdu.writeto(filename)

def _saveJpeg(filename, array2d):
    scipy.misc.imsave(filename, array2d)

def createFlat(imageFileList):
    flatStack = scipy.ndimage.imread(imageFileList[0])
    for imageFile in imageFileList[1:]:
        flatStack += scipy.ndimage.imread(imageFile)
    return flatStack / len(imageFileList)

class DetectedFiber(object):
    def __init__(self, centroidDict):
        self.centroidList = [centroidDict]
        # self.centroids = [pyGuideCentroid]
        # self.imageFiles = [imageFileName]

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
    def rad(self):
        return numpy.average([cent["rad"] for cent in self.centroidList], axis=0, weights=self.counts)

    def belongs2me(self, centroidDict, minSep=MINSEP):
        # if center moves by more than 0.25 pixels
        # doesn't belong
        dist = numpy.linalg.norm(numpy.subtract(centroidDict["xyCtr"], self.xyCtr))
        # if dist < 3:
        #     logging.info("belongs to", dist, self.imageFiles)
        return dist < minSep
        # logging.info("dist!", dist, self.imageFiles)
        # return dist<(self.rad/2.)

    def add2me(self, centroidDict):
        self.centroidList.append(centroidDict)

    def detectedIn(self, imageFileName):
        return imageFileName in self.imageFiles


def multiprocessImage(imageFileList, callFunc, block=False):
    # may want to try map_async
    p = Pool(5)
    if block:
        output = p.map(processImage, imageFileList)
        callFunc(output)
        return None
    else:
        return p.map_async(processImage, imageFileList, callback=callFunc)

def sortDetections(brightestCentroidList, plot=False, minCounts=MINCOUNTS, minSep=MINSEP):
    """Reorganize detection list into groups
    of detections (1 group per fiber)
    """
    detectedFibers = []
    for brightestCentroid in brightestCentroidList:
        isNewDetection = None
        crashMe = False
        if brightestCentroid["counts"] is not None and brightestCentroid["counts"] > MINCOUNTS:
            isNewDetection = True
            # search through every previous detection
            for prevDetection in detectedFibers:
                if prevDetection.belongs2me(brightestCentroid, minSep):
                    if isNewDetection == False:
                        crashMe = True
                        logging.info("bad bad, crash me!")
                    # logging.info("previous detection!!!", brightestCentroid.xyCtr)
                    isNewDetection = False
                    # logging.info("previous detection", os.path.split(imageFile)[-1], prevDetection.imageFiles)
                    prevDetection.add2me(brightestCentroid)
                    # break # assign to the first that works???
                # was this is a new detection?
            if isNewDetection:
                # logging.info('new detection:', os.path.split(imageFile)[-1], brightestCentroid.counts, brightestCentroid.xyCtr)
                detectedFibers.append(DetectedFiber(brightestCentroid))

        if plot:
            imageFile = brightestCentroid["imageFile"]
            color = "r" if isNewDetection else "b"
            fig = plt.figure(figsize=(10,10));plt.imshow(scipy.ndimage.imread(imageFile), vmin=0, vmax=10)#plt.show(block=False)
            plt.scatter(0, 0, s=80, facecolors='none', edgecolors='b')
            if brightestCentroid["xyCtr"] is not None:
               x,y = brightestCentroid["xyCtr"]
               plt.scatter(x, y, s=80, facecolors='none', edgecolors=color)
            frameNumber = int(imageFile.split("img")[-1].split(".")[0])
            zfilled = "%i"%frameNumber
            zfilled = zfilled.zfill(5)
            dd = os.path.split(imageFile)[0]
            nfn = os.path.join(dd, "pyguide%s.png"%zfilled)
            fig.savefig(nfn); plt.close(fig)    # close the figure
        # if crashMe:
        #     raise RuntimeError("Non-unique detection!!!!")
    return detectedFibers

def batchMultiprocess(imageFileDirectory, flatImg, imgBaseName="img", imgExtension="bmp"):
    """! Process all images in given directory
    """
    # detectedFiberList = DetectedFiberList(flatImg)
    imageFiles = glob.glob(os.path.join(imageFileDirectory, "*."+imgExtension))
    nImageFiles = len(imageFiles)
    imageFilesSorted = [os.path.join(imageFileDirectory, "%s%i.%s"%(imgBaseName, num, imgExtension)) for num in range(1,nImageFiles)]
    # warning image files are not sorted as expected, even after explicitly sorting
    # eg 999.jpg > 2000.jpg.  this is bad because image order matters very much
    # furthermore rather than
    # note image files are expected to be 1.jpg, 2.jpg, 3.jpg, ..., 354.jpg...
    # while loop seems weird, but whatever
    frameNumber = frameStartNum
    tstart = time.time()
    brightestCentroidList = multiprocessImage(imageFilesSorted)
    detectedFiberList = sortDetections(brightestCentroidList)
    totaltime = time.time() - tstart
    logging.info("total time", totaltime)
    logging.info(nImageFiles/(totaltime), "frames per second processed")
    logging.info("sorting detections")
    return detectedFiberList

def convToFits(imageFileDirectory, flatImg, frameStartNum, frameEndNum=None, imgBaseName="img", imgExtension="bmp"):
    saveImage()

if __name__ == "__main__":
    def callMeWhenDone(brightestCentroidList):
        detectedFiberList = sortDetections(brightestCentroidList, plot=True)
        logging.info("Done, found ", len(detectedFiberList), "fibers")

    imgDir = "/home/lcomapper/Documents/Camera_test/test061"
    imgBase = "img"
    flatImg = None
    # build the image list (in the correct order, glob doesn't do it correctly)
    unsortedImgs = glob.glob(os.path.join(imgDir, "img*.bmp"))
    frameNum = 1
    sortedImgs = []
    while True:
        imgName = os.path.join(imgDir, "img%i.bmp"%frameNum)
        if imgName in unsortedImgs:
            sortedImgs.append(imgName)
            frameNum += 1
        else:
            break
    logging.info("going to process ", len(sortedImgs), " images")
    block = multiprocessImage(sortedImgs, callMeWhenDone, block=True)

    # flatImgList = [os.path.join(imgDir, "%s%i.bmp"%(imgBase, ii)) for ii in range(nImg-7,nImg)]
    # flatImg = createFlat(flatImgList)
    # frameStartNum = 1
    # frameEndNum = None
    # detectedFiberList = batchProcess(imgDir, flatImg, frameStartNum, frameEndNum, imgBase)
    # detectedFiberList = batchMultiprocess(imgDir, flatImg)
    # remove those not detected in more than 1 frame
    # detectedFiberList = [detectedFiber for detectedFiber in detectedFiberList.detectedFibers if len(detectedFiber.imageFiles)>1]
    # logging.info("Done, found ", len(detectedFiberList.detectedFibers), "fibers")


