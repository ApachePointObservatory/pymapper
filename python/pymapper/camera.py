""" camera control
"""
from __future__ import division, absolute_import

import os
import time
import glob
import subprocess
import collections
import pickle

# import numpy
import scipy.ndimage

from twisted.internet import reactor
from twisted.internet.task import LoopingCall

from pymapper.imgProcess import processImage, multiprocessImage #DetectedFiberList

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
        # self.fileNumSeen = 1
        self.fileNumProc = 1
        # self.unprocessedImages = collections.deque([])
        # self.timeStamps = []
        # self.imgLoadTimes = []
        # self.watchLoop = LoopingCall(self.watchDirectory)
        self.acquiring = False
        self.process = None
        self.imageDir = imageDir
        assert os.path.exists(imageDir), "%s doesn't exist, create it"%imageDir
        # clean up any existing files
        assert len(self.getAllImgFiles())==0, "%s is not empty!"%imageDir
        self.acquisionCB = None
        self.procImgCall = None
        self.centroidList = []

    def doneProcessingCallback(self, callFunc):
        # to be called when all image processing is done,
        # the centroids should have been saved as a picklefile
        self.procImgCall = callFunc

    @property
    def currFile(self):
        return self.getNthFile(self.fileNumProc)

    def getAllImgFiles(self):
        return glob.glob(os.path.join(self.imageDir, "*.%s"%EXTENSION))

    def getNthFile(self, fileNum):
        filename = "%s%i.%s"%(BASENAME, fileNum, EXTENSION)
        return os.path.join(self.imageDir, filename)

    def getRemainingFileList(self):
        nFiles = len(self.getAllImgFiles())
        return [self.getNthFile(fileNum) for fileNum in range(self.fileNumProc, nFiles+1)]

    def beginAcquisition(self, callFunc=None):
        """call callFunc when acquision has started
        """
        if callFunc is not None:
            print("setting acquisiont cb", callFunc)
            self.acquisitionCB = callFunc
        self.acquiring = True
        # self.watchDirectory()
        # self.watchLoop.start(0.) # call repeatedly, non-blocking
        # self.watchDirectory()
        self.process = subprocess.Popen(EXE, cwd=self.imageDir)
        self.processImageLoop()


    def stopAcquisition(self):
        print("Stopping Camera Acquision")
        self.process.kill()
        # self.watchLoop.stop()
        self.acquiring = False
        # calculate frames per second
        # diff = numpy.mean(numpy.diff(self.timeStamps))
        # std = numpy.std(diff)
        # total_secs = self.timeStamps[-1] - self.timeStamps[0]
        # print("fps", len(self.timeStamps)/total_secs)
        # total_secs = self.imgLoadTimes[-1] - self.imgLoadTimes[0]
        # print("fps loaded", len(self.imgLoadTimes)/total_secs)

    # def watchDirectory(self):
    #     currFile = self.getCurrFile()
    #     if os.path.exists(currFile):
    #         print("see file", currFile)
    #         # when the first file is seen, call
    #         # the acquisition callback
    #         self.timeStamps.append(time.time())
    #         # self.unprocessedImages.append(currFile)
    #         if self.fileNumSeen == 1 and self.acquisitionCB is not None:
    #             print("acquisition started")
    #             reactor.callLater(0., self.acquisitionCB)
    #             # self.processImageLoop() # begin processing images
    #             # self.acquisitionCB()


    #         self.fileNumSeen += 1
    #         # free eventloop (paranoia?)
    #         # reactor.callLater(0, self.procImgCall, currFile)
    #     if self.acquiring:
    #         reactor.callLater(0, self.watchDirectory)

    def multiprocessDone(self, remainingCentroidList):
        print("All frames processed!")
        print("pickling centroid list")
        self.centroidList.extend(remainingCentroidList)
        pickleCentroids(self.centroidList, self.imageDir)
        self.procImgCall()

    def processImageLoop(self):
        # called recursively until acquisition is done
        # once acquisition is done switch to multiprocessing
        # the remaining images!

        if os.path.exists(self.currFile):
            if self.fileNumProc == 1 and self.acquisitionCB is not None:
                print("acquisition started")
                reactor.callLater(0., self.acquisitionCB)
            if self.acquiring:
                # continue processing one image at a time
                # recursively here
                self.centroidList.append(processImage(self.currFile))
                self.fileNumProc += 1
                reactor.callLater(0, self.processImageLoop)
                return
            else:
                # done acquiring, move to multiprocessing!
                # blocks?!
                print("beginning multiprocessing")
                self.multiProcessWait = multiprocessImage(self.getRemainingFileList(), self.multiprocessDone)
                return
        else:
            # file wasn't found
            # are we still acquiring?
            if self.acquiring:
                reactor.callLater(0, self.processImageLoop)
                return
            else:
                raise RuntimeError("processImageLoop in weird state")
        # import pdb; pdb.set_trace()

        # try:
        #     fileName = self.getNthFile(self.fileNumProc)
        #     self.procImgCall(fileName)
        #     self.fileNumProc += 1
        # except IOError as e:
        #     print("failed to find a process %s"%fileName)
        #     if not self.acquiring:
        #         # camera is no longer acquiring, must be done
        #         print("camera done! and all frames processed!")
        #         return
        #     # print("e", e)
        # reactor.callLater(0, self.processImageLoop)
        # # if self.unprocessedImages:
        # #     print("%i images queued"%len(self.unprocessedImages))
        # #     self.procImgCall(self.unprocessedImages.popleft())
        # # if self.acquiring or self.unprocessedImages:
        # #     # call me again
        # #     reactor.callLater(0.01, self.processImageLoop)
        # # else:
        # #     print("all images processed")

        # # decide


    def loadImage(self, filename):
        self.imageBuffer.append(scipy.ndimage.imread(filename, flatten=True))
        self.imgLoadTimes.append(time.time())


if __name__ == "__main__":
    camera = Camera()
    camera.start()
    reactor.callLater(3*60., camera.stop)
    reactor.run()