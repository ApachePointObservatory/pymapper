""" Run mapping process
"""
from __future__ import division, absolute_import

import argparse

from math import floor
import os
import glob
import sys
import logging

from twisted.internet import reactor

from sdss.utilities.astrodatetime import datetime

from .camera import Camera, sortDetections, IMGBASENAME, IMGEXTENSION, getScanParams, pickleDetectionList
# from .imgProcess import DetectedFiberList
from .motor import MotorController
from .fiberAssign import SlitheadSolver, FocalSurfaceSolver

homedir = os.path.expanduser("~")
baseDir = os.path.join(homedir, "Documents/Camera_test")
baseName = "test"
fiberslitposfile = os.path.join(os.getenv("PYMAPPER_DIR"), "etc", "fiberpos.dat")


"""
todo: add re-detect/re-solve options
"""

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".

    from: http://stackoverflow.com/questions/3041986/python-command-line-yes-no-input
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

def determineScanNumber(plateID, mjddir):
    # replace plPlugMapP-XXX with plPlugMapM-XXX
    # previousDirectory, filename = os.path.split(self.plPlugMap.filename)
    # determine any existing scans
    globStr = os.path.join(mjddir, "plPlugMapM-%i-*.par"%plateID)
    nExisting = glob.glob(globStr)
    # number this scan accordingly
    return len(nExisting) + 1

def pathPlugMapP(plateID):
    plateZfill = ("%i"%plateID).zfill(6)
    #replace 10s, 1s place with XX
    plateSubDir = plateZfill[:-2] + "XX"
    fileName = "plPlugMapP-%i.par"%plateID
    platelistdir = os.environ["PLATELIST_DIR"]
    return os.path.join(platelistdir, "plates", plateSubDir, plateZfill, fileName)


def getExistingImgs(scanDir):
    imgs = glob.glob(os.path.join(scanDir, "%s*.%s"%(IMGBASENAME, IMGEXTENSION)))
    return imgs

def configureLogging(scanDir, overwrite=True):
    """if overwrite is true, remove any existing logfile,
    else append to any existing one
    """
    # configure logging
    logfile = os.path.join(scanDir, "scan.log")

    if os.path.exists(logfile) and overwrite:
        os.remove(logfile)

    root = logging.getLogger()
    root.setLevel(logging.DEBUG)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # ch.setFormatter(formatter)
    root.addHandler(ch)
    return logfile


def reprocess(args):
    # get all relative information
    # from existing log file
    if args.scanDir is None:
        raise RuntimeError("Must specify --scanDir with --reprocess")
    baseDir = os.path.abspath(args.rootDir)
    scanDir = os.path.join(baseDir, args.scanDir)
    if not os.path.exists(scanDir):
        raise RuntimeError("Scan directory doesn't exist: %s"%scanDir)
    # verify that images exist
    imgs = getExistingImgs(scanDir)
    if not imgs:
        raise RuntimeError("Scan directory doesn't contain existing imgs, cannot --reprocess!")
    logfile = configureLogging(scanDir, overwrite=False)
    logging.info("reprocessing images in %s"%scanDir)
    # determine previous scan params
    scanParams = getScanParams(logfile)
    startPos = scanParams["start"]
    endPos = scanParams["end"]
    scanSpeed = scanParams["speed"]
    # create directory to hold camera images
    # note all previous images will be removed if image dir is not empty
    camera = Camera(scanDir, startPos, scanSpeed)
    def solvePlate():
        # load the (previously pickled centroid list)
        logging.info("sorting detections. makePlots=%s"%str(args.makePlots))
        detectedFiberList = sortDetections(camera.centroidList, plot=args.makePlots)
        # pickle and save the detection list
        pickleDetectionList(detectedFiberList, scanDir)
        logging.info("scan finished.")
        logging.info("found %i fibers"%len(detectedFiberList))
        nImages = len(camera.getAllImgFiles)
        scanTime = abs(endPos - startPos)/float(scanSpeed)
        fps = nImages / scanTime
        logging.info("%i images taken, FPS: %.4f"%(nImages, fps))
        # not yet ready to solve plate:
        logging.info("correlating fiber positions on slit")
        shs = SlitheadSolver(detectedFiberList)
        shs.getOffsetAndScale()
        shs.matchDetections()
        plugMapPath = pathPlugMapP(args.plateID)
        logging.info("plugmap path", plugMapPath)
        assert os.path.exists(plugMapPath)
        fss = FocalSurfaceSolver(detectedFiberList, plugMapPath)
    camera.doneProcessingCallback(solvePlate)
    camera.reprocessImages()


def runScan(args):
    """Move motor, take images, etc
    """
    if args.plateID is None:
        raise RuntimeError("Must specify --plateID")
    baseDir = os.path.abspath(args.rootDir)
        # baseDir doesn't exist
    # verify base directory is ok
    if not os.path.exists(baseDir):
        question = "Base directory %s does not exist, create it?"%baseDir
        output = query_yes_no(question, default="yes")
        if output:
            try:
                os.makedirs(baseDir)
            except:
                logging.info("failed to create base directory %s (permissions?)"%baseDir)
                sys.exit()
        else:
            # exit program
            logging.info("bye!"); sys.exit()

    if args.scanDir is None:
        # no scan dir, try to provide the next one
        # try to determine the next number in the
        # directory sequence
        imgs = sorted(glob.glob(os.path.join(baseDir, "%s*"%baseName)))
        imgInts = []
        for img in imgs:
            try:
                imgInts.append(int(img[-3:]))
            except:
                pass
        if imgInts:
            nextTestInt = max(imgInts)+1
            scanDir = "%s%s"%(baseName, ("%i"%nextTestInt).zfill(3))
        else:
            # no existing test dirs, create the first one?
            scanDir = "test001"
        question = "No scan directory provided, create new one: %s?"%scanDir
        output = query_yes_no(question, default="yes")
        if not output:
            logging.info("bye!"); sys.exit()
    else:
        scanDir = args.scanDir
    scanDir = os.path.join(baseDir, scanDir)
    if os.path.exists(scanDir):
        # scan dir already exists, check if images have been written here yet?
        imgs = getExistingImgs(scanDir)
        if imgs:
            # images exist ask user what to do
            question = "%i pre-existing images in %s, REMOVE ALL?"%(len(imgs), scanDir)
            output = query_yes_no(question, default="yes")
            if output:
                try:
                    for img in imgs:
                        os.remove(img)
                except:
                    logging.info("failed to remove existing images in %s (permissions?)"%scanDir)
                    sys.exit()
    else:
        #create the new scan dir
        try:
            os.makedirs(scanDir)
        except:
            logging.info("failed to create scan directory in %s (permissions?)"%scanDir)
            sys.exit()

    logging.info("scanDir: %s"%scanDir)
    logging.info("plate ID: %i"%args.plateID)
    logging.info("motor start pos (mm): %.2f"%args.startPos)
    logging.info("motor end pos (mm): %.2f"%args.endPos)
    logging.info("motor scan speed (mm/sec): %.2f"%args.scanSpeed)


    # create directory to hold camera images
    # note all previous images will be removed if image dir is not empty
    camera = Camera(scanDir, args.startPos, args.scanSpeed)

    # setup object that finds and holds detections
    # detectedFiberList = DetectedFiberList()

    # construct motor, nothing happens until connect is called
    motorController = MotorController(
        startPos = args.startPos,
        endPos = args.endPos,
        scanSpeed = args.scanSpeed,
        hostname = "10.3.252.13",
        port = 15000,
        )

    # set up callback chains for mapping process

    def stopCamera():
        logging.info("stopCamera")
        # motor is done scanning, kill camera
        camera.stopAcquisition()

    def moveMotor():
        # camera is acquiring begin moving motor/laser
        logging.info("moveMotor")
        motorController.scan(callFunc=stopCamera)

    def startCamera():
        # motor is ready to move, begin snapping pics
        logging.info("startCamera")
        camera.beginAcquisition(callFunc=moveMotor)

    def solvePlate():
        # load the (previously pickled centroid list)
        logging.info("sorting detections. makePlots=%s"%str(args.makePlots))
        detectedFiberList = sortDetections(camera.centroidList, plot=args.makePlots)
        # pickle and save the detection list
        pickleDetectionList(detectedFiberList, scanDir)
        logging.info("scan finished.")
        logging.info("found %i fibers"%len(detectedFiberList))
        nImages = len(glob.glob(os.path.join(scanDir, "*.bmp")))
        scanTime = abs(args.endPos - args.startPos)/float(args.scanSpeed)
        fps = nImages / scanTime
        logging.info("%i images taken, FPS: %.4f"%(nImages, fps))
        shs = SlitheadSolver(detectedFiberList)
        shs.getOffsetAndScale()
        shs.matchDetections()
        # not yet ready to solve plate:

        plugMapPath = pathPlugMapP(args.plateID)
        logging.info("plugmap path", plugMapPath)
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

def main(argv=None):
    global baseDir
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run the plate mapper."
        )
    parser.add_argument("--plateID", type=int, required=False, help="Plate ID")
    parser.add_argument("--scanDir", required=False, help="""Directory in which to put scan.
                            If not provided, one will be automatically determined"""
                            )
    parser.add_argument("--rootDir", required=False,
        default=baseDir, help="Root directory, scanDir will be created here.")
    parser.add_argument("--startPos", required=False, type=float, default=24,
        help="begin of scan motor position (mm).")
    parser.add_argument("--endPos", required=False, type=float, default=134,
        help="end of scan motor position (mm).")
    parser.add_argument("--scanSpeed", required=False, type=float, default=0.6,
        help="speed at which motor scans (mm/sec).")
    parser.add_argument("--makePlots", action="store_true", default=False, help="if present create png plots with circled dectections, takes much longer." )
    parser.add_argument("--solve", action="store_true", default=False, help="if present, solve plate matching fibers to holes.")
    parser.add_argument("--reprocess", action="store_true", default=False, help="if present, reprocess the raw images.")
    args = parser.parse_args()

    if args.reprocess:
        reprocess(args)
    else:
        runScan(args)




if __name__ == "__main__":
    main()
    #try:
    #    main()
    #except Exception as e:
    #    reactor.stop()
    #    raise(e)
       # sys.exit()




"""
    MJD = floor(datetime.now().mjd)
    logging.info("MJD: %i"%MJD)
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
    logging.info("scanNumber %i"%scanNum)
"""
