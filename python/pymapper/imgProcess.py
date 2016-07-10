"""Routines for processing images from the lco mapper
"""
from __future__ import division, absolute_import
import os

import numpy
import scipy.ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


from .camera import getFrameNumFromName


THRESH = 50

def applyThreshold(array2d, thresh):
    # which pixels are greater than sigmaDetect sigma from the mean
    pixInd = numpy.argwhere(array2d < thresh)
    array2d[pixInd] = 0
    return array2d

def to2d(array1d, imshape):
    return numpy.reshape(array1d, imshape)


def createFlat(imageFileList):
    flatStack = scipy.ndimage.imread(imageFileList[0])
    for imageFile in imageFileList[1:]:
        flatStack += scipy.ndimage.imread(imageFile)
    return flatStack / len(imageFileList)



# def batchMultiprocess(imageFileDirectory, flatImg, imgBaseName="img", imgExtension="bmp"):
#     """! Process all images in given directory
#     """
#     # detectedFiberList = DetectedFiberList(flatImg)
#     imageFilesSorted = getSortedImages(imageFileDirectory, imgBaseName, imgExtension)
#     frameNumber = frameStartNum
#     tstart = time.time()
#     brightestCentroidList = multiprocessImage(imageFilesSorted)
#     detectedFiberList = sortDetections(brightestCentroidList)
#     totaltime = time.time() - tstart
#     print("total time", totaltime)
#     print(nImageFiles/(totaltime), "frames per second processed")
#     print("sorting detections")
#     return detectedFiberList



# if __name__ == "__main__":
    # def callMeWhenDone(brightestCentroidList):
    #     detectedFiberList = sortDetections(brightestCentroidList, plot=True)
    #     print("Done, found ", len(detectedFiberList), "fibers")

    # imgDir = "/home/lcomapper/Documents/Camera_test/test061"
    # imgBase = "img"
    # flatImg = None
    # # build the image list (in the correct order, glob doesn't do it correctly)
    # unsortedImgs = glob.glob(os.path.join(imgDir, "img*.bmp"))
    # frameNum = 1
    # sortedImgs = []
    # while True:
    #     imgName = os.path.join(imgDir, "img%i.bmp"%frameNum)
    #     if imgName in unsortedImgs:
    #         sortedImgs.append(imgName)
    #         frameNum += 1
    #     else:
    #         break
    # print("going to process ", len(sortedImgs), " images")
    # block = multiprocessImage(sortedImgs, callMeWhenDone, block=True)

    # flatImgList = [os.path.join(imgDir, "%s%i.bmp"%(imgBase, ii)) for ii in range(nImg-7,nImg)]
    # flatImg = createFlat(flatImgList)
    # frameStartNum = 1
    # frameEndNum = None
    # detectedFiberList = batchProcess(imgDir, flatImg, frameStartNum, frameEndNum, imgBase)
    # detectedFiberList = batchMultiprocess(imgDir, flatImg)
    # remove those not detected in more than 1 frame
    # detectedFiberList = [detectedFiber for detectedFiber in detectedFiberList.detectedFibers if len(detectedFiber.imageFiles)>1]
    # print("Done, found ", len(detectedFiberList.detectedFibers), "fibers")


