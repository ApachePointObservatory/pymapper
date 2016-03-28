""" camera control
"""
from __future__ import division, absolute_import

import os
import time
import glob
import subprocess
import collections

# import numpy
import scipy.ndimage

from twisted.internet import reactor
from twisted.internet.task import LoopingCall

from pymapper.imgProcess import DetectedFiberList

# how to ensure that the image program is killed?

# EXE = "testFakeWrite.py"
# BASENAME = "pyTestFile"
# EXTENSION = "foo"

EXE = "AsynchronousGrabWrite"
BASENAME = "img"
EXTENSION = "bmp"

# self.imageDir = os.path.join(os.environ["PYMAPPER_DIR"], "tests")

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
        for f in glob.glob(os.path.join(self.imageDir, "*.%s"%EXTENSION)):
            os.remove(f)
        # self.imageBuffer = []
        self.acquisionCB = None
        self.procImgCall = None
        self.detectedFiberList = DetectedFiberList()


    def doneProcessingCallback(self, callFunc):
        # to be called when all image processing is done,
        # receives exported detectdFiberList
        self.procImgCall = callFunc

    @property
    def currFile(self):
        return self.getNthFile(self.fileNumProc)

    def getNthFile(self, fileNum):
        filename = "%s%i.%s"%(BASENAME, fileNum, EXTENSION)
        return os.path.join(self.imageDir, filename)

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

    def processImageLoop(self):
        # called recursively until acquisition is done
        # and all images processed
        # if an unprocessed image is on the queue, process it
        # print("unprocessed images", self.unprocessedImages)
        if os.path.exists(self.currFile):
            if self.fileNumProc == 1 and self.acquisitionCB is not None:
                print("acquisition started")
                reactor.callLater(0., self.acquisitionCB)
            self.detectedFiberList.processImage(self.currFile)
            self.fileNumProc += 1
            reactor.callLater(0, self.processImageLoop)
        else:
            # file wasn't found
            # are we still acquiring?
            if self.acquiring:
                reactor.callLater(0, self.processImageLoop)
            else:
                # not acquiring
                print("camera done! and all frames processed!")
                print("found %i fibers"%len(self.detectedFiberList.detectedFibers))
                print("calling matching routines")
                self.procImgCall(self.detectedFiberList.export())
                # import pdb; pdb.set_trace()
                return

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