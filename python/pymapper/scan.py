
from __future__ import division, absolute_import
# import datetime
import glob
import os

import numpy
import scipy.ndimage

# many settings taken from idlmapper/src/evilscan.c

MAX_SCANPIX = 1500000
THRESH = 100       # signal above the mean that we trigger on. Meaningless guess for testing.


# other globals that should come from elsewhere!!!
MOTORSPEED = 400
FPS = 30 # frames per second
MAXMEAN = 50 + 5

class Scan(object):
    def __init__(self, filename):
        """! A Scan

        @param[in] filename, string
        """
        if os.path.exists(filename):
            raise RuntimeError("file already exists, specify another", filename)
        self.filename
        self.writeFileHeader()
        self.frameNum = 0
        self.thresh = None
        self.motor = None

    def writeFileHeader(self):
        """!Initialize the scan file and write header.
        """

        print("\nWriting output file %s\n"%self.filename);
        with open(self.filename, "w") as f:

            # header info, I don't care about this now but
            # put it in eventually

            # f.write("EVILSCAN\n")
            # # f.write("fscanVersion %s\n\n"%scanVersion)
            # f.write("pluggers     %s\n"%pluggers)
            # f.write("plateId      %d\n"%platenum)
            # f.write("fscanMJD     %d\n"%scanMJD)
            # f.write("fscanId      %d\n"%scanId)
            # f.write("fscanDate    %s"%datetime.datetime.now().isoformat()) # includes CR
            # f.write("fscanFile    %s\n\n", filename)
            # f.write("fscanMode    %s\n", scan_mode_long)
            # f.write("fscanSpeed   %d\n", mspeed)
            # f.write("fscanRows    %d\n", rows)
            # f.write("fscanCols    %d\n", cols)
            # f.write("fscanBias    %f\n", bias)
            # f.write("motorId1     %d\n", motorID1)
            # f.write("motorId2     %d\n", motorID2)
            # f.write("motorId3     %d\n\n", motorID3)
            # f.write("cartridgeId  %d\n\n", cartid)

            f.write("typedef struct {\n")
            f.write("  int    motor\n")
            f.write("  int    frame\n")
            f.write("  int    motorpos\n")
            f.write("  float  tstamp\n")
            f.write("  int    row\n")
            f.write("  int    col\n")
            f.write("  int    flux\n")
            f.write("} SCANPIX\n\n")

    def processFrame(self, imageFile):
        """! Process a single image

        @param[in] imageFile. String image to be opened)
        """
        imgData = scipy.ndimage.imread(imageFile, flatten=True)
        if self.thresh is None:
            # use this frame (the first to determine the bias level)
            # note we could compute bias level every time without
            # much overhead...
            self.thresh = numpy.mean(imgData) + THRESH
        imgMean = numpy.mean(imgData)
        if imgMean > MAXMEAN:
             raise RuntimeError("TOO MUCH SIGNAL FROM SCATTERED LIGHT ON THE CAMERA!")
        litUpPixels = numpy.nonzero(imgData > self.thresh)
        if len(litUpPixels) > MAX_SCANPIX:
            raise RuntimeError("TOO MUCH SIGNAL FROM SCATTERED LIGHT ON THE CAMERA!")
        timestamp = self.frameNum * FPS
        motorPos = timestamp * MOTORSPEED
        with open(self.filename, "a") as f:
            for litUpPixel in litUpPixels:
                dataLine = "SCANPIX %2d %5d %5d %9.4f %3d %3d %3d\n"%(
                    self.motor,
                    self.frame,
                    motorPos,
                    timestamp,
                    litUpPixel[0],
                    litUpPixel[1],
                    imgData[litUpPixel],
                    )
                f.write(dataLine)
        self.frameNum += 1

    def batchProcess(self, imageFileDirectory, motorID, imgExtension="jpeg"):
        """! Process all images in given directory, assigned to motorID

        @param[in]: imageFileDirectory: directory containing image files (in sortable order)
        @param[in]: motorID: 1, 2, or 3
        """
        imageFiles = glob.glob(os.path.join(imageFileDirectory, "*."+imgExtension))
        self.motor = motorID
        self.frameNum = 0
        for imageFile in imageFiles:
            self.processFrame(imageFileDirectory)


