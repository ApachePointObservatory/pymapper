"""Routines for processing images from the lco mapper
"""
from __future__ import division, absolute_import
import glob
import os
import collections
from operator import attrgetter
import time
import pickle

import numpy
import scipy.ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from astropy.io import fits

import PyGuide

THRESH = 50
MINCOUNTS = 50

def applyThreshold(array2d, thresh):
    # which pixels are greater than sigmaDetect sigma from the mean
    pixInd = numpy.argwhere(array2d < thresh)
    array2d[pixInd] = 0
    return array2d

def pyGuideFind(array2d, medianFilter=False):
    if medianFilter:
        array2d = scipy.ndimage.median_filter(array2d, size=2)
    ccdInfo = PyGuide.CCDInfo(bias=50, readNoise=10, ccdGain=1)
    findStars = PyGuide.findStars(array2d, None, None, ccdInfo)
    print("found ", len(findStars[0]), "fibers")
    plt.figure(figsize=(13,13))
    plt.imshow(array2d)
    for centroid in findStars[0]:
        print("centroid counts:", centroid.counts)
        x,y = centroid.xyCtr
        plt.scatter(x, y, s=80, facecolors='none', edgecolors='r')
    plt.show()
    return findStars

def to2d(array1d, imshape):
    return numpy.reshape(array1d, imshape)

def saveImage(filename, array2d):
    if os.path.exists(filename):
        print("removing previous", filename)
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
    def __init__(self, pyGuideCentroid, imageFrameName):
        self.centroids = [pyGuideCentroid]
        self.imageFrames = [imageFrameName]

    @property
    def xyCtr(self):
        return self.centroids[-1].xyCtr

    @property
    def rad(self):
        return self.centroids[-1].rad

    def belongs2me(self, pyGuideCentroid):
        dist = numpy.linalg.norm(numpy.subtract(pyGuideCentroid.xyCtr, self.xyCtr))
        return dist<(self.rad/2.)

    def add2me(self, pyGuideCentroid, imgFrameName):
        self.centroids.append(pyGuideCentroid)
        self.imageFrames.append(imgFrameName)

    def detectedIn(self, imageFrameName):
        return imageFrameName in self.imageFrames

class DetectedFiberList(object):
    def __init__(self, flatImg):
        self.ccdInfo = PyGuide.CCDInfo(bias=50, readNoise=10, ccdGain=1)
        # detected fibers is keyed by imageName
        self.detectedFibers = []
        self.flatImg = flatImg
        self.prevFrames = collections.deque([], maxlen=4)

    def pickle(self, fileName):
        detectedFibers = []
        for detectedFiber in self.detectedFibers:
            detectedFibers.append(
                dict((
                    ("imageFrames", detectedFiber.imageFrames),
                    ("counts", [centroid.counts for centroid in detectedFiber.centroids]),
                    ("xyCtrs", [centroid.xyCtr for centroid in detectedFiber.centroids]),
                ))
            )

        # save detections to picked file
        output = open(fileName, "wb")
        pickle.dump(detectedFibers, output)
        output.close()

    @property
    def lastFrameDetections(self):
        """Return all detected fibers found in the previous image frame
        """
        previousFibers = []
        for prevImg in self.prevFrames:
            for fiber in self.detectedFibers[::-1]:
                if fiber.detectedIn(prevImg):
                    previousFibers.append(fiber)
                else:
                    break
        return previousFibers

    def processImage(self, imageFile, frameNumber):
        """! Process a single image

        @param[in] imageFile. String

        """
        # read in the image data and apply the flat
        imgData = scipy.ndimage.imread(imageFile) #/ self.flatImg
        # imgData = scipy.ndimage.gaussian_filter(imgData, 1, mode="nearest")
        #plt.figure();plt.imshow(imgData);plt.show()
        pyGuideCentroids = PyGuide.findStars(imgData, None, None, self.ccdInfo)[0]
        brightestCentroid = pyGuideCentroids[0]
        # fig = plt.figure(figsize=(10,10));plt.imshow(imgData, vmin=0, vmax=10)#plt.show(block=False)
        # plt.scatter(0, 0, s=80, facecolors='none', edgecolors='b')
        # for centroid in pyGuideCentroids:
        #     x,y = centroid.xyCtr
        #     plt.scatter(x, y, s=80, facecolors='none', edgecolors='r')
        # zfilled = "%i"%frameNumber
        # zfilled = zfilled.zfill(5)
        # dd = os.path.split(imageFile)[0]
        # nfn = os.path.join(dd, "pyguide%s.png"%zfilled)
        # fig.savefig(nfn); plt.close(fig)    # close the figure


        # print("max, mean value: ", numpy.max(imgData), numpy.mean(imgData))
        # toss all finds with counts below the threshold:
        # pyGuideCentroids = [pyGuideFind for pyGuideFind in pyGuideCentroids if pyGuideFind.counts > MINCOUNTS]
        # determine if any of these detections were present in the previous frame,
        # only look to the previous image (fibers may only be detected in contiguous images)
        # so don't look further back than one image.
        # if so, apply them to the correct (previous) detection
        # if they were not previously detected, create a new detection
        # print("found ", len(pyGuideCentroids), "fibers ")
        # prevDetections = self.lastFrameDetectionsf
        # for ind, pyGuideFind in enumerate(pyGuideCentroids):
            # print("found", pyGuideFind.xyCtr, pyGuideFind.counts, pyGuideFind.rad)
            # is this a new dectection or was it found already in the previous image?
        isNewDetection = True
        for prevDetection in self.lastFrameDetections: #self.detectedFibers:
            if prevDetection.belongs2me(brightestCentroid):
                if isNewDetection == False:
                    raise RuntimeError("This shouldn't ever happen")
                # print("previous detection!!!", brightestCentroid.xyCtr)
                isNewDetection = False
                prevDetection.add2me(brightestCentroid, imageFile)
                # break # assign to the first that works???
            # was this is a new detection?
        if isNewDetection:
            print('new detection: img', imageFile)
            self.detectedFibers.append(DetectedFiber(brightestCentroid, imageFile))
        self.prevFrames.append(imageFile)

def batchProcess(imageFileDirectory, flatImg, frameStartNum, frameEndNum=None, imgBaseName="img", imgExtension="bmp"):
    """! Process all images in given directory
    """
    detectedFiberList = DetectedFiberList(flatImg)
    imageFiles = glob.glob(os.path.join(imageFileDirectory, "*."+imgExtension))
    # warning image files are not sorted as expected, even after explicitly sorting
    # eg 999.jpg > 2000.jpg.  this is bad because image order matters very much
    # furthermore rather than
    # note image files are expected to be 1.jpg, 2.jpg, 3.jpg, ..., 354.jpg...
    # while loop seems weird, but whatever
    frameNumber = frameStartNum
    tstart = time.time()
    while True:
        nextImageFile = "%s%i.%s"%(imgBaseName, frameNumber, imgExtension)
        imageFilePath = os.path.join(imageFileDirectory, nextImageFile)
        if imageFilePath not in imageFiles:
            print(imageFilePath, "not in imageFiles", len(imageFiles))
            break
        # print("processing: ", imageFilePath)
        detectedFiberList.processImage(imageFilePath, frameNumber)
        if frameEndNum is not None and frameEndNum == frameNumber:
            break
        frameNumber += 1
    print((frameNumber-frameStartNum)/(time.time()-tstart), "frames per second processed")
    return detectedFiberList

def convToFits(imageFileDirectory, flatImg, frameStartNum, frameEndNum=None, imgBaseName="img", imgExtension="bmp"):
    saveImage()

if __name__ == "__main__":
    imgDir = "/home/lcomapper/Desktop/mappertestingMarch/run5"
    imgBase = "img"
    flatImgList = [os.path.join(imgDir, "%s%i.bmp"%(imgBase, ii)) for ii in range(1,7)]
    flatImg = createFlat(flatImgList)
    frameStartNum = 10
    frameEndNum = None
    detectedFiberList = batchProcess(imgDir, flatImg, frameStartNum, frameEndNum, imgBase)
    # remove those not detected in more than 1 frame
    # detectedFiberList = [detectedFiber for detectedFiber in detectedFiberList.detectedFibers if len(detectedFiber.imageFrames)>1]
    print("Done, found ", len(detectedFiberList.detectedFibers), "fibers")

    # save detections to picked file
    pickleFile = os.path.join(imgDir, "pickledDetections.pkl")
    detectedFiberList.pickle(pickleFile)


    # import pdb; pdb.set_trace(
    # for fiber in detectedFiberList.detectedFibers: #sorted(detectedFiberList.detectedFibers, key=attrgetter("xyCtr")):
    #     print (fiber.xyCtr, [os.path.split(x)[-1] for x in fiber.imageFrames])


    # imgNumber = 83
    # imgData = scipy.ndimage.imread(os.path.join(imgDir, "%s%i.bmp"%(imgBase, imgNumber)))
    # # pyGuideFind(imgData)
    # saveImage(os.path.join(imgDir, "%s%i.fits"%(imgBase, imgNumber)), imgData)
    # # now check the flat
    # flatImg = imgData/flatImg
    # print("flat applied!")
    # pyGuideFind(flatImg)
    # saveImage(os.path.join(imgDir, "%s%i_proc.fits"%(imgBase, imgNumber)), flatImg)

