from __future__ import division, absolute_import

import itertools
import copy

import scipy.ndimage

# from sdss.utilities.yanny import yanny

class FScanFrame(object):
    def __init__(self, frame, flatImg, motor, motorpos, tstamp, row, col, flux):
        """
        Class representing a single fscan frame
        @param[in] frame: int, frame number
        @param[in] flatImg: 2D array. A flat image from the CCD.
        @param[in] motor: int, motor ID (1-3)
        @param[in] tstamp: float, timestamp for frame
        @param[in] row: int, row of pixel
        @param[in] col: int, col of pixel
        @param[in] flux: int, flux of pixel
        """
        self.frame = frame
        self.flatImg = flatImg
        self.motor = motor
        self.motorpos = motorpos
        self.tstamp = tstamp
        self.rows = [row]
        self.cols = [col]
        self.fluxes = [flux]

    def addPixel(self, row, col, flux):
        """Add a single pixel and it's flux

        @param[in] row: int, row of pixel
        @param[in] col: int, col of pixel
        @param[in] flux: int, flux of pixel
        """
        self.rows.append(row)
        self.cols.append(col)
        self.fluxes.append(flux)

    def getImg(self):
        """Return a 2D image of this frame
        Add fluxes from specfic pixels to the flat image
        """
        imgData = copy.deepcopy(self.flatImg)
        for row, col, flux in itertools.izip(self.rows, self.cols, self.fluxes):
            imgData[col, row] += flux
        return imgData

class ScanMovie(object):
    def __init__(self, fscanFile, flatImgFile):
        """Class for representing an evil mapper fscan as a pymapper movie

        @param[in]  fscanFile: path to fscan file
        @param[in]  flatImgFile: path to a flatImgFile
        """
        self.fscanFile = fscanFile
        self.flatImg = scipy.ndimage.imread(flatImgFile)
        self.frames = self.framesFromFile(fscanFile)

    def parseLine(self, line):
        """Parse a single line beginning with SCANPIX
        """
        spix, motor, frame, motorpos, tstamp, row, col, flux = line.split()
        return (int(motor), int(frame), int(motorpos),
            float(tstamp), int(row), int(col), int(flux))

    def framesFromFile(self, fscanFile):
        """Parse the fscan file return an array of frame objects
        """
        frames = []
        with open(fscanFile) as f:
            filelines = f.readlines()
        currentFrame = None
        for line in filelines:
            if not line.startswith("SCANPIX"):
                continue
            motor, frame, motorpos, tstamp, row, col, flux = self.parseLine(line)
            if currentFrame is None:
                # first frame
                currentFrame = FScanFrame(frame, self.flatImg, motor, motorpos, tstamp, row, col, flux)
            elif currentFrame.frame != frame:
                # new frame
                frames.append(currentFrame)
                currentFrame = FScanFrame(frame, self.flatImg, motor, motorpos, tstamp, row, col, flux)
            else:
                # append to current frame
                currentFrame.addPixel(row, col, flux)
        # add final frame
        frames.append(currentFrame)
        return frames



if __name__ == "__main__":
    fscanFile = "data/8637/fiberScan-8637-57534-01.par"
    flatFile = "data/flat.bmp"
    sm = ScanMovie(fscanFile, flatFile)
    import pdb; pdb.set_trace()