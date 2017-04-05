#!/usr/bin/env python
import os
import glob
import shutil

from pymapper import runMapper

def parsePlPlug(plPlugPath):
    outDict = {
        "plateId" : None,
        "cartridgeId" : None,
        "fscanMJD" : None,
        "fscanId" : None,
    }
    with open(plPlugPath, "r") as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip
        for key in outDict.keys():
            if line.startswith(key):
                value = int(line.split())
                outDict[key] = value
    assert not None in outDict.values()
    return outDict


fromDir = "/data/rawmapper/57843/plate8785/fscan9"
walkDirs = os.walk(fromDir)
for path, dirs, files in walkDirs:
    if "fscan" in path:
        print("working with ", path)
        plPlug = glob.glob(path + "/plPlug*.par")
        if not len(plPlug) == 1:
            print("no unique plPlugmatch", plPlug, "skipping")
            continue
        try:
            scanPars = parsePlPlug(plPlug[0])
        except:
            print("failed to parse", plPlug[0])
            continue
        print("moving old scan data to saved directory")
        savedPath = os.path.join(path, "prefiberswap")
        os.makedirs(savedPath)
        tocopy = glob.glob(path + "/*.png") + glob.glob(path + "/*.par") + glob.glob("/*.pkl")
        for copyme in tocopy:
            base, filename = os.path.split()
            newfile = os.path.join(savedPath, filename)
            shutil.copyfile(copyme, newfile)



if __name__ == "__main__":
    runMapper.resolve()
