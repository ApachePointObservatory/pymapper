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

testDir = os.path.join(os.environ["PYMAPPER_DIR"], "tests")

class Camera(object):
    def __init__(self):
        self.fileNum = 1
        self.timeStamps = []
        self.imgLoadTimes = []
        self.watchLoop = LoopingCall(self.watchDirectory)
        self.process = None
        # clean up any existing files
        for f in glob.glob(os.path.join(testDir, "*.%s"%EXTENSION)):
            os.remove(f)
        self.imageBuffer = []

    def getCurrFile(self):
        return self.getNthFile(self.fileNum)

    def getNthFile(self, fileNum):
        filename = "%s%i.%s"%(BASENAME, fileNum, EXTENSION)
        return os.path.join(testDir, filename)

    def start(self):
        self.watchLoop.start(0.) # call repeatedly, non-blocking
        self.process = subprocess.Popen(EXE, cwd=testDir)

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
            self.timeStamps.append(time.time())
            self.fileNum += 1
            # free eventloop (paranoia?)
            reactor.callLater(0, self.loadImage, currFile)

    def loadImage(self, filename):
        self.imageBuffer.append(scipy.ndimage.imread(filename, flatten=True))
        self.imgLoadTimes.append(time.time())


if __name__ == "__main__":
    camera = Camera()
    camera.start()
    reactor.callLater(3*60., camera.stop)
    reactor.run()