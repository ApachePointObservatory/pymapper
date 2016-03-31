""" Run mapping process
"""
from __future__ import division, absolute_import

from math import floor
import os
import glob

from twisted.internet import reactor

from sdss.utilities.astrodatetime import datetime

from pymapper.imgProcess import sortDetections
from pymapper.camera import Camera#, unpickleCentroids
# from pymapper.imgProcess import DetectedFiberList
from pymapper.motor import MotorController
from pymapper.fiberAssign import FocalSurfaceSolver

homedir = os.path.expanduser("~")
scandir = os.path.join(homedir, "scan")
platelistdir = os.environ["PLATELIST_DIR"]



def determineScanNumber(plateID, mjddir):
    # replace plPlugMapP-XXX with plPlugMapM-XXX
    # previousDirectory, filename = os.path.split(self.plPlugMap.filename)
    # determine any existing scans
    globStr = os.path.join(mjddir, "plPlugMapM-%i-*.par"%plateID)
    print("globStr", globStr)
    nExisting = glob.glob(globStr)
    # number this scan accordingly
    return len(nExisting) + 1

def pathPlugMapP(plateID):
    plateZfill = ("%i"%plateID).zfill(6)
    #replace 10s, 1s place with XX
    plateSubDir = plateZfill[:-2] + "XX"
    fileName = "plPlugMapP-%i.par"%plateID
    return os.path.join(platelistdir, "plates", plateSubDir, plateZfill, fileName)


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
    plateDir = os.path.join(mjddir, "%i"%plateID)
    if not os.path.exists(plateDir):
        os.makedirs(plateDir)
    # create a scan number dir
    scanNum = 1
    while True:
        imageDir = os.path.join(plateDir, "rawImage-%i"%scanNum)
        if not os.path.exists(imageDir):
            os.makedirs(imageDir)
            break
        scanNum += 1
    print("scanNumber %i"%scanNum)
    # create directory to hold camera images
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

    def solvePlate():
        # load the (previously pickled centroid list)
        detectedFiberList = sortDetections(camera.centroidList)
        plugMapPath = pathPlugMapP(plateID)
        print("plugmap path", plugMapPath)
        assert os.path.exists(plugMapPath)
        fss = FocalSurfaceSolver(detectedFiberList, plugMapPath)

    camera.doneProcessingCallback(solvePlate)

    motorController.addReadyCallback(startCamera)
    # hand fiber detection routine to the camera, so it is called
    # when any new frame is available
    # camera.addProcImageCall(detectedFiberList.processImage)
    # connecting to the motor starts it all off
    motorController.connect()
    reactor.run()

 #   import pickle
 #   pkl = open(os.path.join(imageDir, "detectionList.pkl"), "rb")
 #   detectedFiberList = pickle.load(pkl)
 #   pkl.close()
 #   solvePlate(detectedFiberList)




