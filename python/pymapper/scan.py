"""More or less copied from idlmapper/src/evilscan.c
"""
from __future__ import division, absolute_import
from sdss.utilities.astrodatetime import datetime
import glob
import os

import numpy
import scipy.ndimage

__all__ = ["Scan"]

MAX_SCANPIX = 1500000
THRESH = 50 #100       # signal above the mean that we trigger on. Meaningless guess for testing.


# other globals that should come from elsewhere , eg configuration file!!!
MOTORSPEED = 400
FPS = 30 # frames per second
MAXMEAN = 50 + 5
# hacks to make evilmap4 work
MOTORID1 = 46
MOTORID2 = 31
MOTORID3 = 14
# CARTID = 1
CAM_HEIGHT = 960
CAM_WIDTH = 960
SCANID = 1

class Scan(object):
    def __init__(self, plate, filename, motorSpeed=MOTORSPEED, fps=FPS, mjd=None):
        """! A Scan

        @param[in] plate, int
        @param[in] filename, string
        @param[in] motorSpeed, float
        @param[in] fps, frames per second, float
        @param[in] mjd, mdj or None
        """
        self.plate = plate
        self.motorSpeed = motorSpeed
        self.fps = fps
        if os.path.exists(filename):
            raise RuntimeError("file %s already exists, specify another, or remove it"%filename)
        self.filename = filename
        self.thresh = None
        self.bias = None
        self.motor = None
        self.now = datetime.now()
        self.mjd = mjd if mjd is not None else self.now.mjd
        self.dataLines = [] # will hold SCANPIX lines to be written to a file

    def writeOutputFile(self):
        """!Write the scan file
        """

        print("\nWriting output file %s\n"%self.filename)
        with open(self.filename, "w") as fileObj:
            # header info, I don't care about this now but
            # put it in eventually

            # fileObj.write("EVILSCAN\n")
            # # fileObj.write("fscanVersion %s\n\n"%scanVersion)
            # fileObj.write("pluggers     %s\n"%plateuggers)
            fileObj.write("plateId      %d\n"%self.plate)
            fileObj.write("fscanMJD     %d\n"%self.mjd) # needed for evilmap4
            fileObj.write("fscanId      %d\n"%SCANID) # needed for evilmap4
            fileObj.write("fscanDate    %s\n"%self.now.isoformat()) # includes CR
            fileObj.write("fscanFile    %s\n\n"%self.filename)
            # fileObj.write("fscanMode    %s\n", scan_mode_long)
            fileObj.write("fscanSpeed   %d\n"%self.motorSpeed)
            fileObj.write("fscanRows    %d\n"%CAM_WIDTH) # needed for evilmap4
            fileObj.write("fscanCols    %d\n"%CAM_HEIGHT) # needed for evilmap4
            fileObj.write("fscanBias    %f\n"%self.bias) # needed for evilmap4
            fileObj.write("motorId1     %d\n"%MOTORID1) # needed for evilmap4
            fileObj.write("motorId2     %d\n"%MOTORID2) # needed for evilmap4
            fileObj.write("motorId3     %d\n\n"%MOTORID3) # needed for evilmap4
            # fileObj.write("cartridgeId  %d\n\n", cartid)

            fileObj.write("typedef struct {\n")
            fileObj.write("  int    motor\n")
            fileObj.write("  int    frame\n")
            fileObj.write("  int    motorpos\n")
            fileObj.write("  float  tstamp\n")
            fileObj.write("  int    row\n")
            fileObj.write("  int    col\n")
            fileObj.write("  int    flux\n")
            fileObj.write("} SCANPIX\n\n")

            # write collected lines of data
            for line in self.dataLines:
                fileObj.write(line)


    def processFrame(self, imageFile, frameNumber, motorID):
        """! Process a single image

        @param[in] imageFile. String
        @param[in] frameNumber. Int

        """
        imgData = scipy.ndimage.imread(imageFile, flatten=True)
        if self.thresh is None and self.bias is None:
            # use the 5th frame to determine the bias
            # use this frame (the first to determine the bias level)
            # note we could compute bias level every time without
            # much overhead...
            self.bias = numpy.mean(imgData)
            self.thresh = self.bias + THRESH
        imgMean = numpy.mean(imgData)
        if imgMean > MAXMEAN:
             raise RuntimeError("TOO MUCH SIGNAL FROM SCATTERED LIGHT ON THE CAMERA!")
        litUpPixels = numpy.argwhere(imgData > self.thresh)
        if len(litUpPixels) > MAX_SCANPIX:
            raise RuntimeError("TOO MUCH SIGNAL FROM SCATTERED LIGHT ON THE CAMERA!")
        timestamp = frameNumber / self.fps
        motorPos = timestamp * self.motorSpeed
        for litUpPixel in litUpPixels:
            x,y = litUpPixel
            dataLine = "SCANPIX %2d %5d %5d %9.4f %3d %3d %3d\n"%(
                motorID,
                frameNumber,
                motorPos,
                timestamp,
                x,
                y,
                imgData[x,y],
                )
            self.dataLines.append(dataLine)

    def batchProcess(self, imageFileDirectory, motorID, imgExtension="jpg"):
        """! Process all images in imageFileDirectory

        @param[in]: imageFileDirectory: directory containing image files (in sortable order)
        @param[in]: motorID: 1, 2, or 3
        @param[in]: imgExtension, used for globbing.
        """
        imageFiles = glob.glob(os.path.join(imageFileDirectory, "*."+imgExtension))
        # warning image files are not sorted as expected, even after explicitly sorting
        # eg 999.jpg > 2000.jpg.  this is bad because image order matters very much
        # note image files are expected to be 1.jpg, 2.jpg, 3.jpg, ..., 354.jpg...
        # while loop seems weird, but whatever
        frameNumber = 1
        while True:
            nextImageFile = "%i.%s"%(frameNumber, imgExtension)
            imageFilePath = os.path.join(imageFileDirectory, nextImageFile)
            if imageFilePath not in imageFiles:
                print("all images done, %i frames processed"%frameNumber)
                break
            if frameNumber%100==0:
                # print progress
                print("%.1f percent done, motor %i"%(100*frameNumber/float(len(imageFiles)), motorID))
            self.processFrame(imageFilePath, frameNumber, motorID)
            frameNumber += 1

