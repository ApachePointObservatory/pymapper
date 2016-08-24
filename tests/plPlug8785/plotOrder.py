from __future__ import division, absolute_import
import itertools

import numpy
import matplotlib.pyplot as plt

from sdss.utilities.yanny import yanny

# unmapped = "plPlug8785P2.par"
unmapped = "plPlugMapP-8785.par"
mapped = "test169plug.par"
mapped2 = "../data/test141/plPlugMapM-8785-55555-01-test.par"
overlay2pl = "fiberBlocksAPOGEE_SOUTH.unordered.par"
# map1 = "plPlugMapM_1.par"
# map2 = "plPlugMapM_2.par"
# map3 = "plPlugMapM_3.par"
# map4 = "plPlugMapM_4.par"
# map5 = "plPlugMapM_5.par"
# map6 = "plPlugMapM_6.par"
# map7 = "plPlugMapM_7.par"
mappedList = ["plPlugMapM_%i.par"%ii for ii in range(1,8)]

#plate 8785
# unmappedPlPlug = yanny(filename=unmapped, np=True)
# objectInds = numpy.argwhere(unmappedPlPlug["PLUGMAPOBJ"]["holeType"]=="OBJECT")
# xPos = unmappedPlPlug["PLUGMAPOBJ"]["xFocal"][objectInds].flatten()
# yPos = unmappedPlPlug["PLUGMAPOBJ"]["yFocal"][objectInds].flatten()

plLookup = yanny(filename=overlay2pl, np=True)
fiberIDs = plLookup["TIFIBERBLOCK_APOGEE_SOUTH"]["fiberid"] # blanton says this is mapping between plplug (index) and overlay (number in array)
# blanton's map isn't right use my photoshop by eye matches

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


# lookupX = plLookup["TIFIBERBLOCK_APOGEE_SOUTH"]["fibercenx"]
# lookupY = plLookup["TIFIBERBLOCK_APOGEE_SOUTH"]["fiberceny"]

# plt.figure()
# for ii, (fiberID, x, y) in enumerate(itertools.izip(fiberIDs, lookupX, lookupY)):
#     plt.text(x, y, str(fiberID), horizontalalignment="center", fontsize=10, color="green")
#     plt.text(x, y+0.1, str(ii), horizontalalignment="center", fontsize=10, color="blue")
# plt.xlim([min(lookupX), max(lookupX)])
# plt.ylim([min(lookupY), max(lookupY)])
# plt.show()

def rmsErrMappedVsUnmapped():
    plPlugMap = yanny(filename=unmapped, np=True)
    objectInds = numpy.argwhere(plPlugMap["PLUGMAPOBJ"]["holeType"]=="OBJECT")
    xPosUnmapped = numpy.asarray(plPlugMap["PLUGMAPOBJ"]["xFocal"][objectInds].flatten())
    yPosUnmapped = numpy.asarray(plPlugMap["PLUGMAPOBJ"]["yFocal"][objectInds].flatten())
    plPlugMap = yanny(filename=mapped, np=True)
    objectInds = numpy.argwhere(plPlugMap["PLUGMAPOBJ"]["holeType"]=="OBJECT")
    xPosMapped = numpy.asarray(plPlugMap["PLUGMAPOBJ"]["xFocal"][objectInds].flatten())
    yPosMapped = numpy.asarray(plPlugMap["PLUGMAPOBJ"]["yFocal"][objectInds].flatten())
    rms = numpy.sum(numpy.sqrt((xPosUnmapped-xPosMapped)**2+(yPosUnmapped-yPosMapped)**2))
    print("RMS: %.5f"%rms)

def plotUnmappedIndexXY():
    fig = plt.figure(figsize=(50,50))
    plLookup = yanny(filename=overlay2pl, np=True)
    fiberIDs = plLookup["TIFIBERBLOCK_APOGEE_SOUTH"]["fiberid"] # blanton says this is mapping between plplug (index) and overlay (number in array)
    plPlugMap = yanny(filename=unmapped, np=True)
    objectInds = numpy.argwhere(plPlugMap["PLUGMAPOBJ"]["holeType"]=="OBJECT")
    xPos = numpy.asarray(plPlugMap["PLUGMAPOBJ"]["xFocal"][objectInds].flatten())
    yPos = numpy.asarray(plPlugMap["PLUGMAPOBJ"]["yFocal"][objectInds].flatten())
    # invert x
    xPos = xPos*-1
    for ii, (x,y,blantonInd) in enumerate(itertools.izip(xPos,yPos, fiberIDs)):
        plt.text(x,y,str(ii+1), horizontalalignment='center', fontsize=22, color="black")
        plt.text(x,y+3,str(blantonInd), horizontalalignment='center', fontsize=22, color="red")
    plt.xlim([min(xPos), max(xPos)])
    plt.ylim([min(yPos), max(yPos)])
    fig.savefig("blantonVsPlugMap-8785.png"); plt.close(fig)



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

    # expectedOrder = []
    # mappedOrder = []
    # for ii, ind in enumerate(objectInds):
    #     expectedOrder.append(ii)
    #     mappedOrder.append(plPlugMap["PLUGMAPOBJ"]["fiberId"][ind])
    # return expectedOrder, mappedOrder

def reorderInBlock(mappedOrder):
    # slice into blocks of 30, and remap backwards

    blocks = numpy.asarray(mappedOrder).reshape((10,30))
    flipped = blocks[:,::-1]
    flattened = flipped.flatten()
    return flattened

def reorderBlocks(mappedOrder):
    # slice into blocks of 30, and remap backwards
    blocks = numpy.asarray(mappedOrder).reshape((10,30))
    flipped = blocks[::-1,:]
    flattened = flipped.flatten()
    return flattened

def reorderAll(mappedOrder):
    blocks = numpy.asarray(mappedOrder).reshape((10,30))
    flipped = blocks[::-1,::-1]
    flattened = flipped.flatten()
    return flattened

def reorderCustom(mappedOrder):
    blocks = numpy.asarray(mappedOrder).reshape((10,30))
    # flip order and direction of first 3 blocks
    blocks[:3,:] = blocks[:3,:][::-1,::-1]
    # flip order and directon of next 2 blocks
    blocks[3:5,:] = blocks[3:5,:][::-1,::-1]
    # flip order and director of next 3 blocks
    blocks[5:8,:] = blocks[5:8,:][::-1,::-1]
    # flip order and direction of last set of blocks
    blocks[8:,:] = blocks[8:,:][::-1,::-1]
    flattened = blocks.flatten()
    return flattened

def _playWithOrder(mappedOrder):
    blocks = numpy.asarray(mappedOrder).reshape((10,30))
    blocksOut = numpy.array(blocks, copy=True)*0-10
    blocksOut[0,:] = blocks[8,:]
    blocksOut[1,:] = blocks[9,:]
    blocksOut[2,:] = blocks[5,:]
    blocksOut[3,:] = blocks[6,:]
    blocksOut[4,:] = blocks[7,:]
    blocksOut[5,:] = blocks[3,:]
    blocksOut[6,:] = blocks[4,:]
    blocksOut[7,:] = blocks[0,:]
    blocksOut[8,:] = blocks[1,:]
    blocksOut[9,:] = blocks[2,:]
    plt.figure()
    for ii, block in enumerate(blocksOut):
        if ii == 0:
            print(sorted(block))
        expectedNums = numpy.arange(ii*30, ii*30+30) + 1
        plt.plot(expectedNums, block, 'o-')
    # blocksOut[5,:] = blocks[2,:]
    plt.show()

lookupDict = numpy.argsort(numpy.asarray(range(25,31) + range(17,25) + range(9,17) + range(1,9))-1)
lookupDict = numpy.argsort(numpy.asarray(range(8,0, -1) + range(16,8,-1) + range(24,16,-1) + range(30,24,-1))-1)

def playWithOrder(mappedOrder, expectedOrder):
    blocks = numpy.asarray(mappedOrder).reshape((10,30))
    order = numpy.asarray(expectedOrder).reshape((10,30))
    # blocksOut = numpy.array(blocks, copy=True)*0-10
    # blocksOut[0,:] = blocks[8,:]
    # blocksOut[1,:] = blocks[9,:]
    # blocksOut[2,:] = blocks[5,:]
    # blocksOut[3,:] = blocks[6,:]
    # blocksOut[4,:] = blocks[7,:]
    # blocksOut[5,:] = blocks[3,:]
    # blocksOut[6,:] = blocks[4,:]
    # blocksOut[7,:] = blocks[0,:]
    # blocksOut[8,:] = blocks[1,:]
    # blocksOut[9,:] = blocks[2,:]
    plt.figure()
    for ii, block in enumerate(blocks):
        if ii == 0:
            print(sorted(block))
        # for jj,b in itertools.izip(lookupDict,block):
        expectedNums = numpy.arange(ii*30, ii*30+30) + 1
        # expectedNums = expectedNums[lookupDict]
            # plt.plot(jj-1, b, 'o')
        # expectedNums = order[ii,:]
        plt.plot(expectedNums, block, 'o-')
    # blocksOut[5,:] = blocks[2,:]
    plt.show()

# plotUnmappedIndexXY()
# rmsErrMappedVsUnmapped()
# plt.figure()
eo, mo = getOrder(mapped)
# # mo = mo[numpy.argsort(eo)]
plt.plot(eo, mo, 'ob', alpha=0.5)
# eo, mo = getOrder(mapped2)
# plt.plot(eo, mo, 'o-r', alpha=0.5)

# playWithOrder(mo, fiberIDs)
# expededOrder = fiberIDs
# plt.plot(eo, mo, 'o-r', alpha=0.5)

# plt.plot(eo, reorderBottomUp(mo), 'o-k')

# plt.figure()
# plt.hist(mo, 1200)
# plt.xlim([-2,302])
# plt.ylim([0, 1.1])
# plt.plot(eo, reorderAll(mo), 'o-b')
plt.show()
