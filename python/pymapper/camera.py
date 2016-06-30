""" camera control
"""
from __future__ import division, absolute_import

import os
import glob
import subprocess
import pickle
import time
import logging

from twisted.internet import reactor

from pymapper.imgProcess import multiprocessImage
# how to ensure that the image program is killed?

# EXE = "testFakeWrite.py"
# BASENAME = "pyTestFile"
# EXTENSION = "foo"

EXE = "AsynchronousGrabWrite"
BASENAME = "img"
EXTENSION = "bmp"

# self.imageDir = os.path.join(os.environ["PYMAPPER_DIR"], "tests")

def pickleCentroids(centroidList, imageDir):
    # save detections to picked file
    fileName = os.path.join(imageDir, "centroidList.pkl")
    output = open(fileName, "wb")
    pickle.dump(centroidList, output)
    output.close()

def unpickleCentroids(imageDir):
   pkl = open(os.path.join(imageDir, "centroidList.pkl"), "rb")
   centroidList = pickle.load(pkl)
   pkl.close()
   return centroidList

class Camera(object):
    def __init__(self, imageDir):
        self.acquiring = False
        self.process = None
        self.imageDir = imageDir
        assert os.path.exists(imageDir), "%s doesn't exist, create it"%imageDir
        # whine if this directory already has images
        assert len(self.getAllImgFiles())==0, "%s is not empty!"%imageDir
        self.acquisionCB = None
        self.procImgCall = None
        self.centroidList = []

    def doneProcessingCallback(self, callFunc):
        # to be called when all image processing is done,
        # the centroids should have been saved as a picklefile
        self.procImgCall = callFunc

    def getAllImgFiles(self):
        # glob doesn't order correctly in
        return glob.glob(os.path.join(self.imageDir, "*.%s"%EXTENSION))

    def getNthFile(self, fileNum):
        filename = "%s%i.%s"%(BASENAME, fileNum, EXTENSION)
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
        if callFunc is not None:
            print("setting acquision cb", callFunc)
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
        print("processed %.2f frames per second"% (len(self.centroidList)/(time.time()-self.processingStart)))
        pickleCentroids(self.centroidList, self.imageDir)
        if self.procImgCall is not None:
            self.procImgCall()

    def waitForFirstImage(self):
        # loop here until the first image is seen.  Once it is
        # fire the acquisition callback and begin processing
        if os.path.exists(self.getNthFile(1)):
            print("acquisition started")
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

    # def multiprocessNext(self, centroidList):
    #     print("multiprocessNext")
    #     self.centroidList.extend(centroidList)
    #     self.multiprocessImageLoop()
        # reactor.callLater(0., self.multiprocessImageLoop)

    def multiprocessImageLoop(self, centroidList=None):
        # called recursively until all images are (multi!) processed
        print("multiprocessImageLoop")
        if centroidList:
            print("adding %i centroids"%len(centroidList))
            self.centroidList.extend(centroidList)
        unprocessedFileList = self.getUnprocessedFileList()
        # don't process more than 20 images at a time
        unprocessedFileList = unprocessedFileList[:50]
        if unprocessedFileList:
            print("processing images %s to %s"%tuple([os.path.split(_img)[-1] for _img in [unprocessedFileList[0], unprocessedFileList[-1]]]))
            nonBlock = multiprocessImage(unprocessedFileList, self.multiprocessImageLoop, block=False)
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



if __name__ == "__main__":
    camera = Camera()
    camera.start()
    reactor.callLater(3*60., camera.stop)
    reactor.run()