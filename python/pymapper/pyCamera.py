from __future__ import absolute_import, division, print_function
import os
import xml.etree.ElementTree as ET
import re
import numpy
from astropy.io import fits
import time
import pymba
import shutil

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

class Timer(object):
    def __init__(self):
        self.tstart = numpy.nan
    def elapsed(self):
        return time.time() - self.tstart

timer = Timer()

configFile = os.path.join(os.getenv("PYMAPPER_DIR"), "etc/mapperCamConfig.xml")

cameraConfig = VimbaConfig(configFile)
cameraSettings = cameraConfig.allSettings

EXPTIME = cameraSettings["ExposureTimeAbs"]
GAIN = cameraSettings["GainRaw"]
if cameraSettings["PixelFormat"] == "Mono8":
    numpydtype = numpy.uint8
else:
    numpydtype = numpy.uint16
frameNum = 0

# imgSaveDir = os.path.join(os.path.expanduser("~"), "tmpMapImg")
imgSaveDir = None

# if os.path.exists(imgSaveDir):
#     # delete it to ensure it is empty!
#     shutil.rmtree(imgSaveDir)
# os.makedirs(imgSaveDir)

def frameCB(frame):
    global frameNum
    imgData = numpy.ndarray(buffer = frame.getBufferByteData(),
                           dtype = numpydtype,
                           shape = (frame.height,
                                    frame.width)
                            )
    frameNum += 1
    if frameNum % 100 == 0:
        print("fps: %.4f, imgNum: %i"%(frameNum/timer.elapsed(), frameNum))
    strNum = ("%i"%frameNum).zfill(6)
    filename = os.path.join(imgSaveDir, "img%s.fits"%strNum)
    hdu = fits.PrimaryHDU(imgData)
    hdulist = fits.HDUList([hdu])
    prihdr = hdulist[0].header
    prihdr["tstamp"] = time.time(), "UNIX time of exposure"
    prihdr["exptime"] = EXPTIME, "EXPTIME micro seconds"
    prihdr["gain"] = GAIN, "GAIN value in decibels"
    hdulist.writeto(filename)
    hdulist.close()
    frame.queueFrameCapture(frameCB)
    # del hdu
    # del hdulist
    # del prihdr

def loadConfig(camera):
    """Explicitly set all configuation specified in the config obj
    """
    for key, val in cameraSettings.iteritems():
        print("setting: %s = %s"%(key, str(val)))
        setattr(camera, key, val)

vimba = pymba.Vimba()
vimba.startup()
system = vimba.getSystem()
system.runFeatureCommand("GeVDiscoveryAllOnce")
camera = vimba.getCamera(cameraConfig.cameraID)
camera.openCamera()
loadConfig(camera)
frames = [
    camera.getFrame(),
    camera.getFrame(),
    camera.getFrame()
    ]
for frame in frames:
    frame.announceFrame()
    frame.queueFrameCapture(frameCB)
camera.startCapture()

def startCapture(imgSaveDir):
    """capture shortcut, don't save
    """
    global imgSaveDir
    imgSaveDir = imgSaveDir
    timer.tstart = time.time()
    camera.runFeatureCommand("AcquisitionStart")

def stopCapture():
    timer.tstart = numpy.nan
    camera.runFeatureCommand("AcquisitionStop")
    camera.endCapture()
    camera.revokeAllFrames()

if __name__ == "__main__":
    imgSaveDir = os.path.join(os.path.expanduser("~"), "tmpMapImg")
    if os.path.exists(imgSaveDir):
        # delete it to ensure it is empty!
        shutil.rmtree(imgSaveDir)
    os.makedirs(imgSaveDir)
    startCapture(imgSaveDir)
    time.sleep(10)
    stopCapture()

