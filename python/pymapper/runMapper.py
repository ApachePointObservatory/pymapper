""" Run mapping process
"""
from __future__ import division, absolute_import
# import shutil

from functools import partial
import subprocess

import argparse
import time

from math import floor
import os
import glob
import sys

from twisted.internet import reactor

from sdss.utilities.astrodatetime import datetime

from .camera import Camera, sortDetections, IMGBASENAME, IMGEXTENSION, getScanParams, \
    pickleDetectionList, unpickleCentroids
# from .imgProcess import DetectedFiberList
from .motor import MotorController
from .fiberAssign import SlitheadSolver, FocalSurfaceSolver

from . import plt

MOTOR_IP = "10.1.1.26"
MOTOR_PORT = 15000
MJD = floor(datetime.now().mjd + 0.4) # MJD + 0.4 convention chosen by Holtz
BASEDIR = os.path.join(os.path.expanduser("~"), "scan", "%i"%MJD)
if not os.path.exists(BASEDIR):
    os.makedirs(BASEDIR)
# check for existing scan directories

fiberslitposfile = os.path.join(os.getenv("PYMAPPER_DIR"), "etc", "fiberpos.dat")
tstart = None

# import cProfile, pstats, StringIO
# pr = cProfile.Profile()

"""
todo: add re-detect/re-solve options
"""

def loadPlPlugMapM(plPlugMapPath):
    """Loads a plPlugMapM file to the DB.
    largely stolen from sdss_python_module/bin/plPlugMapM !!!
    """

    from sdss.internal.database.connections import LCODatabaseAdminLocalConnection
    from sdss.internal.database.apo.platedb.plPlugMapM import PlPlugMapMFile

    if not os.path.exists(plPlugMapPath):
        raise RuntimeError('file {0} cannot be found'.format(fScan))

    plFile = PlPlugMapMFile(plPlugMapPath, verbose=True)
    plFile.load(replace=True, active=True)

    return

class LogStdOut(object):
    def __init__(self, filename):
        """Simply logs all standard out to a log
        """
        self.terminal = sys.stdout
        self.log = open(filename, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass

class LogStdErr(object):
    def __init__(self, filename):
        """Simply logs all standard out to a log
        """
        self.terminal = sys.stderr
        self.log = open(filename, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass

def configureLogging(scanDir, MJD, plateID, fscanID):
    logfile = os.path.join(scanDir, "scan-%i-%i-fscan%i.log"%(MJD, plateID, fscanID))
    errfile = os.path.join(scanDir, "scan_err-%i-%i-fscan%i.log"%(MJD, plateID, fscanID))
    sys.stdout = LogStdOut(logfile)
    sys.stderr = LogStdErr(errfile)

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

# def configureLogging(scanDir, overwrite=True):
#     """if overwrite is true, remove any existing logfile,
#     else append to any existing one
#     """
#     # configure logging
#     logfile = os.path.join(scanDir, "scan.log")

#     if os.path.exists(logfile) and overwrite:
#         os.remove(logfile)

#     root = logging.getLogger()
#     root.setLevel(logging.DEBUG)

#     ch = logging.StreamHandler(sys.stdout)
#     ch.setLevel(logging.DEBUG)
#     # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
#     # ch.setFormatter(formatter)
#     root.addHandler(ch)
#     return logfile

def _solvePlate(scanDir, plateID, cartID, fscanID, fscanMJD, plot=False, plugMapPath=None, dbLoad=True):
    # global pr
    # pr.disable()
    # s = StringIO.StringIO()
    # sortby = 'cumulative'
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    # print(s.getvalue())
    global tstart
    centroidList = unpickleCentroids(scanDir)
    detectedFiberList = sortDetections(centroidList, plot=plot)
    pickleDetectionList(detectedFiberList, scanDir)
    shs = SlitheadSolver(detectedFiberList, centroidList)
    shs.matchDetections()
    print("missing fibers: ")
    for fiber in shs.missingFibers:
        print("fiber %i"%fiber)
    print("slit match rms: %.8f"%shs.rms)
    shs.plotSolution(scanDir)
    if plugMapPath is None:
        plugMapPath = pathPlugMapP(plateID)
    print("plugmap path", plugMapPath)
    assert os.path.exists(plugMapPath)
    fss = FocalSurfaceSolver(detectedFiberList, plugMapPath, scanDir, cartID, fscanID, fscanMJD)
    print("mapping done: took %.2f seconds"%(time.time()-tstart))
    if shs.missingFibers:
        subprocess.call("gnome-open %s/unplugged.png"%(scanDir), shell = True)
    # load the plPlugMap file as an active plugging in db
    if dbLoad:
        print("Loading plPlugMap in db, and making it active!!!")
        loadPlPlugMapM(fss.plPlugMap.filePath)

        #LCO HACK!!!! un-comment when write access returns
        print("WARNING!!!! NOT COPYING TO /DATA! FIX when perms are set correctly")

        #print("copying plPlugMap file to /data/mapper/<MJD>")
        #basePath, fileName = os.path.split(fss.plPlugMap.filePath)
        #shutil.copy(fss.plPlugMap.filePath, os.path.join("/data/mapper/%i"%MJD, fileName))
        subprocess.call(["cp", fss.plPlugMap.filePath, os.path.join("/data/mapper/%i"%MJD, fileName)])
        print("killing all python processes")
        #print("closing screen log")
        #subprocess.call("exit")
        subprocess.call("killall -9 python", shell=True)
        # plt.show()

# def resolvePlate(args):
#     if args.scanDir is None:
#         raise RuntimeError("Must specify --scanDir with --resolvePlate")
#     baseDir = os.path.abspath(args.rootDir)
#     scanDir = os.path.join(baseDir, args.scanDir)
#     if not os.path.exists(scanDir):
#         raise RuntimeError("Scan directory doesn't exist: %s"%scanDir)
#     _solvePlate(scanDir, args.plateID, plot=args.plotDetections, plugMapPath=args.plPlugMap)


def reprocessImgs(args):
    # global pr
    # pr.enable()
    # get all relative information
    # from existing log file
    if args.scanDir is None:
        raise RuntimeError("Must specify --scanDir with --reprocessImgs")
    baseDir = os.path.abspath(args.rootDir)
    scanDir = os.path.join(baseDir, args.scanDir)
    if not os.path.exists(scanDir):
        raise RuntimeError("Scan directory doesn't exist: %s"%scanDir)
    # verify that images exist
    imgs = getExistingImgs(scanDir)
    if not imgs:
        raise RuntimeError("Scan directory doesn't contain existing imgs, cannot --reprocess!")
    # logfile = configureLogging(scanDir, overwrite=False)
    print("reprocessing images in %s"%scanDir)
    scanParams = getScanParams(os.path.join(scanDir, "scanParams.par"))
    startPos = scanParams["start"]
    endPos = scanParams["end"]
    scanSpeed = scanParams["speed"]
    # create directory to hold camera images
    # note all previous images will be removed if image dir is not empty
    camera = Camera(scanDir, startPos, endPos, scanSpeed)
    solvePlate = partial(_solvePlate, scanDir=scanDir, plateID=args.plateID, plot=args.plotDetections)
    camera.doneProcessingCallback(solvePlate)
    camera.reprocessImages()


def runScan(args):
    """Move motor, take images, etc
    """
    if not args.plateID:
        plateID = int(raw_input("plateID: "))
    else:
        plateID = int(args.plateID)
    if not args.cartID:
        cartID = int(raw_input("cartID: "))
    else:
        cartID = int(args.cartID)
    plateDir = os.path.join(BASEDIR, "plate%i"%plateID)
    if not os.path.exists(plateDir):
        os.makedirs(plateDir)
    # now determine how many scans already exist for this plate
    for fscanID in range(1,10000):
        scanDir = os.path.join(plateDir, "fscan%i"%fscanID)
        if not os.path.exists(scanDir):
            break
    os.makedirs(scanDir)
    # if the MJD directory in /data/mapper doesn't yet exist, make it now
    # plPlugMaps will also be copied there (for utah to grab)

    # LCOHACK!!!!
    print("NOT CREATING /data/MJD dir, no permissions!")
    if not os.path.exists("/data/mapper/%i"%MJD):
        os.makedirs("/data/mapper/%i"%MJD)
    #configureLogging(scanDir, MJD, plateID, fscanID)
    # begin logging screen
    #subprocess.Popen("script %s"%(os.path.join(scanDir, "scan.log")), shell=True)
    print("scanDir: %s"%scanDir)
    print("plate ID: %i"%plateID)
    print("motor start pos (mm): %.2f"%args.startPos)
    print("motor end pos (mm): %.2f"%args.endPos)
    print("motor scan speed (mm/sec): %.2f"%args.scanSpeed)


    # create directory to hold camera images
    # note all previous images will be removed if image dir is not empty
    camera = Camera(scanDir, args.startPos, args.endPos, args.scanSpeed)

    # setup object that finds and holds detections
    # detectedFiberList = DetectedFiberList()

    # construct motor, nothing happens until connect is called
    motorController = MotorController(
        startPos = args.startPos,
        endPos = args.endPos,
        scanSpeed = args.scanSpeed,
        hostname = MOTOR_IP,
        port = MOTOR_PORT,
        )

    # set up callback chains for mapping process

    moveMotor = partial(motorController.scan, callFunc=camera.stopAcquisition)
    beginImgAcquisition = partial(camera.beginAcquisition, callFunc=moveMotor)
    solvePlate = partial(_solvePlate, scanDir, plateID, cartID, fscanID, MJD, plot=args.plotDetections)
    camera.doneProcessingCallback(solvePlate)
    motorController.addReadyCallback(beginImgAcquisition)
    motorController.connect()
    reactor.run()

def recentroid(argv=None):
    """Rerun a map from existing image files
    """
    pass

def redetect(argv=None):
    """Rerun a map from existing centroidList
    """
    pass

def resolve(scanDir=None):
    """Rerun a map from an existing detectionList, finds expected input files (images, centriods, detections, plPlugMapP)
    from the working directory, writes output to the same directory (images, plPlugMapM)
    """
    global tstart
    tstart = time.time()
    if scanDir is None:
        scanDir = os.getcwd()
    plPlugMapFile = glob.glob(os.path.join(scanDir, "plPlugMapP-*.par"))
    if not plPlugMapFile:
        raise RuntimeError("No plPlugMapP file found")
    if not len(plPlugMapFile)==1:
        raise RuntimeError("Found multiple plPlugMap files!")
    plPlugMapFile = plPlugMapFile[0]
    # get the plateID from the plPlug filename
    plateID = int(os.path.split(plPlugMapFile)[-1].strip("plPlugMapP-").strip(".par"))
    _solvePlate(scanDir, plateID, 99, 99, 9999, plot=False, plugMapPath=plPlugMapFile, dbLoad=False)


def main(argv=None):
    global baseDir
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run the plate mapper."
        )
    parser.add_argument("--plateID", type=int, required=False, help="Plate ID")
    parser.add_argument("--cartID", type=int, required=False, help="Cart ID")
    parser.add_argument("--startPos", required=False, type=float, default=134,
        help="begin of scan motor position (mm).")
    parser.add_argument("--endPos", required=False, type=float, default=24,
        help="end of scan motor position (mm).")
    parser.add_argument("--scanSpeed", required=False, type=float, default=1.2,
        help="speed at which motor scans (mm/sec).")
    parser.add_argument("--plotDetections", action="store_true", default=False, help="if present create png plots with circled dectections, takes much longer." )
    # parser.add_argument("--resolvePlate", action="store_true", default=False, help="if present, solve plate matching fibers to holes, from pickled img process output.")
    # parser.add_argument("--reprocessImgs", action="store_true", default=False, help="if present, reprocess the raw images.")
    # parser.add_argument("--plPlugMap", required=False, help="path to plPlugMap file, if not present try to get it from $PLATELIST_DIR")
    args = parser.parse_args()
    global tstart
    tstart = time.time()
    runScan(args)




if __name__ == "__main__":
    main()
