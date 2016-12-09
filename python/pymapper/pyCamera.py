from __future__ import absolute_import, division, print_function
import os
import xml.etree.ElementTree as ET
import re
import numpy
from astropy.io import fits
import time
import pymba
import shutil

# import PyGuide

# # CCDInfo = PyGuide.CCDInfo(bias=50, readNoise=10, ccdGain=1)
# CCDInfo = PyGuide.CCDInfo(bias=1, readNoise=0, ccdGain=1)
# # pyguide mask
# PyGuideMask = numpy.zeros((960,960))
# midPix = 960/2.
# for x in range(960):
#     for y in range(960):
#         pixRad = numpy.sqrt((x-midPix)**2+(y-midPix)**2)
#         if pixRad > midPix:
#             PyGuideMask[x,y]=1

class VimbaConfig(object):
    RW = "R/W"
    Integer = "Integer"
    Boolean = "Boolean"
    Float = "Float"
    Enumeration = "Enumeration"
    def __init__(self, xmlFile):
        """An object representing a configuration state of an AlliedVision GigE camera

        @param[in] xmlFile: configuation file to be loaded (the output of VimbaViewer configuation save)
        """
        self.integerSettings = {}
        self.booleanSettings = {}
        self.enumSettings = {}
        self.floatSettings = {}
        self.fromXML(xmlFile)

    @property
    def allSettings(self):
        return dict(self.integerSettings.items() + self.booleanSettings.items() + self.enumSettings.items() + self.floatSettings.items())

    def fromXML(self, xmlFile):
        """populate settings dicts from xmlFile

        @param[in] xmlFile: configuation file to be loaded (the output of VimbaViewer configuation save)
        """
        # Vimba-written xml files are not parseable without this hack:
        # http://stackoverflow.com/questions/38853644/python-xml-parseerror-junk-after-document-element
        with open(xmlFile) as f:
            xml = f.read()
        tree = ET.fromstring(re.sub(r"(<\?xml[^>]+\?>)", r"\1<root>", xml) + "</root>")
        cameraSettings = tree.find("CameraSettings")
        # save this camera ID
        self.cameraID = cameraSettings.attrib["CameraID"]
        # iterate over all (relevant features in xml file)
        for feature in cameraSettings.find("FeatureGroup").getchildren():
            # only keep read/write values
            if feature.attrib["Access"] != self.RW:
                continue
            value = feature.text
            name = feature.attrib["Name"]
            featureType = feature.attrib["Type"]
            # cast value into expected type
            if featureType == self.Enumeration:
                # leave as a string
                self.enumSettings[name] = str(value)
            elif featureType == self.Integer:
                self.integerSettings[name] = int(value)
            elif featureType == self.Boolean:
                self.booleanSettings[name] = value == "True"
            elif featureType == self.Float:
                self.floatSettings[name] = float(value)
            else:
                print("unknown feature type: %s: %s"%(name, featureType))
                continue


class Globals(object):
    def __init__(self):
        self.configFile = os.path.join(os.getenv("PYMAPPER_DIR"), "etc/mapperCamConfig.xml")
        self.cameraConfig = VimbaConfig(self.configFile)
        self.cameraSettings = self.cameraConfig.allSettings
        self.tstart = -1
        self.imgSaveDir = None
        self.exptime = self.cameraSettings["ExposureTimeAbs"]
        self.gain = self.cameraSettings["GainRaw"]
        if self.cameraSettings["PixelFormat"] == "Mono8":
            self.numpydtype = numpy.uint8
        else:
            self.numpydtype = numpy.uint16
        self.frameNum = 0

    def elapsed(self):
        return time.time() - self.tstart

GLOBALS = Globals()
# if os.path.exists(imgSaveDir):
#     # delete it to ensure it is empty!
#     shutil.rmtree(imgSaveDir)
# os.makedirs(imgSaveDir)

def frameCB(frame):
    imgData = numpy.ndarray(buffer = frame.getBufferByteData(),
                           dtype = GLOBALS.numpydtype,
                           shape = (frame.height,
                                    frame.width)
                            )
    # pyGuideCentroids = PyGuide.findStars(imgData, PyGuideMask, None, CCDInfo)[0]
    GLOBALS.frameNum += 1
    elapsedTime = GLOBALS.elapsed()
    if GLOBALS.frameNum % 100 == 0:
        print("fps: %.4f, imgNum: %i"%(GLOBALS.frameNum/elapsedTime, GLOBALS.frameNum))
    strNum = ("%i"%GLOBALS.frameNum).zfill(6)
    filename = os.path.join(GLOBALS.imgSaveDir, "img%s.fits"%strNum)
    hdu = fits.PrimaryHDU(imgData)
    hdulist = fits.HDUList([hdu])
    prihdr = hdulist[0].header
    prihdr["tstamp"] = elapsedTime, "UNIX time of exposure"
    prihdr["exptime"] = GLOBALS.exptime, "EXPTIME micro seconds"
    prihdr["gain"] = GLOBALS.gain, "GAIN value in decibels"
    prihdr["fnum"] = GLOBALS.frameNum, "GAIN value in decibels"
    hdulist.writeto(filename)
    hdulist.close()
    frame.queueFrameCapture(frameCB)
    # del hdu
    # del hdulist
    # del prihdr

vimba = pymba.Vimba()
vimba.startup()
system = vimba.getSystem()
system.runFeatureCommand("GeVDiscoveryAllOnce")
camera = vimba.getCamera(GLOBALS.cameraConfig.cameraID)
camera.openCamera()
# load config settings
for key, val in GLOBALS.cameraSettings.iteritems():
    print("setting: %s = %s"%(key, str(val)))
    setattr(camera, key, val)
frames = [
    camera.getFrame(),
    camera.getFrame(),
    camera.getFrame(),
    ]

for frame in frames:
    frame.announceFrame()
    frame.queueFrameCapture(frameCB)
camera.startCapture()

def startCapture(imgSaveDir):
    """capture shortcut, don't save
    """
    print("starting image capture")
    GLOBALS.imgSaveDir = imgSaveDir
    GLOBALS.tstart = time.time()
    camera.runFeatureCommand("AcquisitionStart")

def stopCapture():
    print("stopping image capture")
    GLOBALS.tstart = -1
    camera.runFeatureCommand("AcquisitionStop")
    camera.endCapture()
    camera.revokeAllFrames()
    vimba.shutdown()


if __name__ == "__main__":
    imgSaveDir = os.path.join(os.path.expanduser("~"), "tmpMapImg")
    if os.path.exists(imgSaveDir):
        # delete it to ensure it is empty!
        shutil.rmtree(imgSaveDir)
    os.makedirs(imgSaveDir)
    startCapture(imgSaveDir)
    time.sleep(10)
    stopCapture()

