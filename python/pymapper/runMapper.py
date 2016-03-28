""" Run mapping process
"""
from __future__ import division, absolute_import

from math import floor
import os
import glob

from twisted.internet import reactor

from sdss.utilities.astrodatetime import datetime

from pymapper.camera import Camera
# from pymapper.imgProcess import DetectedFiberList
from pymapper.motor import MotorController

homedir = os.path.expanduser("~")
scandir = os.path.join(homedir, "scan")

def determineScanNumber(plateID, mjddir):
    # replace plPlugMapP-XXX with plPlugMapM-XXX
    # previousDirectory, filename = os.path.split(self.plPlugMap.filename)
    # determine any existing scans
    globStr = os.path.join(mjddir, "plPlugMapM-%i-*.par"%plateID)
    print("globStr", globStr)
    nExisting = glob.glob(globStr)
    # number this scan accordingly
    return len(nExisting) + 1

if __name__ == "__main__":

    operator = str(raw_input("Enter Name(s): "))
    cartID = int(raw_input("Enter Cart Number: "))
    plateID = int(raw_input("Enter Plate Number: "))

    MJD = floor(datetime.now().mjd)
    print("MJD: %i"%MJD)
    #create mjd directory
    mjddir = os.path.join(scandir, "%i"%MJD)
    if not os.path.exists(mjddir):
        os.makedirs(mjddir)
    # determine which scan number (dont overwrite any existing)
    scanNumber = determineScanNumber(plateID, mjddir)
    print("scanNumber %i"%scanNumber)
    # create directory to hold camera images
    imageDir = os.path.join(mjddir, "rawImage-%i-%i-%i"%(plateID, MJD, scanNumber))
    if not os.path.exists(imageDir):
        os.makedirs(imageDir)
    # note all previous images will be removed if image dir is not empty
    camera = Camera(imageDir)

    # setup object that finds and holds detections
    # detectedFiberList = DetectedFiberList()

    # construct motor, nothing happens until connect is called
    motorController = MotorController()

    # set up callback chains for mapping process

    def stopCamera():
        print("stopCamera")
        # motor is done scanning, kill camera
        camera.stopAcquisition()

    def moveMotor():
        # camera is acquiring begin moving motor/laser
        print("moveMotor")
        motorController.scan(callFunc=stopCamera)

    def startCamera():
        # motor is ready to move, begin snapping pics
        print("startCamera")
        camera.beginAcquisition(callFunc=moveMotor)

    motorController.addReadyCallback(startCamera)
    # hand fiber detection routine to the camera, so it is called
    # when any new frame is available
    # camera.addProcImageCall(detectedFiberList.processImage)
    # connecting to the motor starts it all off
    motorController.connect()
    reactor.run()



