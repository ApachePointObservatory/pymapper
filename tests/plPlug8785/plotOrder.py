from __future__ import division, absolute_import
import itertools

import numpy
import matplotlib.pyplot as plt

from sdss.utilities.yanny import yanny


unmapped = "plPlugMapP-8785.par"
mapped = "test170plug.par" # output from mapping 8785 determinstically


def getMapFromOverlay():
    # return a list of indices that
    # map the plPlugMap xy positions
    # to the printed overlay xy positions
    with open("overlay2plugMap.txt", "r") as f:
        lines = f.readlines()
    overlayinds = []
    plInds = []
    for line in lines:
        line.strip()
        if line.startswith("#"):
            continue
        overlayInd, plInd = [int(ind) for ind in line.split()]
        if overlayinds:
            assert overlayInd > overlayinds[-1]
        overlayinds.append(overlayInd)
        plInds.append(plInd)

    std = sorted(plInds)
    prev = std[0]
    for ind in std[1:]:
        if ind <= prev:
            print("bad pl ind:", ind)
        prev = ind


    assert len(set(range(1,301)) - set(overlayinds)) == 0
    assert len(set(range(1,301)) - set(plInds)) == 0
    return plInds



def getOrder(plPlugMapFile, invertInts=False):
    # invert ints should be true for an unmapped file
    plPlugMap = yanny(filename=plPlugMapFile, np=True)
    objectInds = numpy.argwhere(plPlugMap["PLUGMAPOBJ"]["holeType"]=="OBJECT")
    xPos = plPlugMap["PLUGMAPOBJ"]["xFocal"][objectInds].flatten()
    yPos = plPlugMap["PLUGMAPOBJ"]["yFocal"][objectInds].flatten()
    radPos = numpy.sqrt(xPos**2+yPos**2)
    plateID = int(plPlugMap["plateId"])
    fiberIDs = plPlugMap["PLUGMAPOBJ"]["fiberId"][objectInds].flatten()
    if invertInts:
        fiberIDs = [-1*fiberID for fiberID in fiberIDs]
    # expectedOrder = numpy.arange(len(fiberIDs))+1
    expectedOrder = numpy.arange(300)[numpy.argsort(getMapFromOverlay())]
    mappedOrder = fiberIDs
    print("eo", expectedOrder[0], expectedOrder[-1])
    return expectedOrder, mappedOrder





eo, mo = getOrder(mapped)
plt.plot(eo, mo, 'ob', alpha=0.5)
plt.show()
