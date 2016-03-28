""" camera control
"""
from __future__ import division, absolute_import

import os
import time
import glob
import subprocess

import numpy
import scipy.ndimage

from twisted.internet import reactor
from twisted.internet.task import LoopingCall

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
        self.fileNum = 1
        self.timeStamps = []
        self.imgLoadTimes = []
        self.watchLoop = LoopingCall(self.watchDirectory)
        self.process = None
        self.imageDir = imageDir
        assert os.path.exists(imageDir), "%s doesn't exist, create it"%imageDir
        # clean up any existing files
        for f in glob.glob(os.path.join(self.imageDir, "*.%s"%EXTENSION)):
            os.remove(f)
        self.imageBuffer = []
        self.acquisionCB = None
        self.procImgCall = None

    def addProcImageCall(self, callFunc):
        # callFunc receives path to file
        self.procImgCall = callFunc

    def getCurrFile(self):
        return self.getNthFile(self.fileNum)

    def getNthFile(self, fileNum):
        filename = "%s%i.%s"%(BASENAME, fileNum, EXTENSION)
        return os.path.join(self.imageDir, filename)

    def beginAcquisition(self, callFunc=None):
        """call callFunc when acquision has started
        """
        if callFunc is not None:
            self.acquisitionCB = callFunc
        self.watchLoop.start(0.) # call repeatedly, non-blocking
        self.process = subprocess.Popen(EXE, cwd=self.imageDir)


    def stop(self):
        self.process.kill()
        self.watchLoop.stop()
        # calculate frames per second
        diff = numpy.mean(numpy.diff(self.timeStamps))
        std = numpy.std(diff)
        total_secs = self.timeStamps[-1] - self.timeStamps[0]
        print("fps", len(self.timeStamps)/total_secs)
        total_secs = self.imgLoadTimes[-1] - self.imgLoadTimes[0]
        print("fps loaded", len(self.imgLoadTimes)/total_secs)

    def watchDirectory(self):
        currFile = self.getCurrFile()
        if os.path.exists(currFile):
            # when the first file is seen, call
            # the acquisition callback
            if self.fileNum == 1 and self.acquisionCB is not None:
                reactor.callLater(0., self.acquisitionCB)
            self.timeStamps.append(time.time())
            self.fileNum += 1
            # free eventloop (paranoia?)
            reactor.callLater(0, self.procImgCall, currFile)

    def loadImage(self, filename):
        self.imageBuffer.append(scipy.ndimage.imread(filename, flatten=True))
        self.imgLoadTimes.append(time.time())


if __name__ == "__main__":
    camera = Camera()
    camera.start()
    reactor.callLater(3*60., camera.stop)
    reactor.run()