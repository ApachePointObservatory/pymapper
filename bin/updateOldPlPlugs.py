#!/usr/bin/env python
import os
import glob
import shutil
import subprocess

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
        line = line.strip()
        for key in outDict.keys():
            if line.startswith(key):
                value = int(line.split(None)[-1])
                outDict[key] = value
    assert not None in outDict.values()
    return outDict


fromDir = "/data/rawmapper/57799/plate9652/fscan1"
walkDirs = os.walk(fromDir)
for path, dirs, files in walkDirs:
    if "fscan" in path and not "prefiberswap" in path:
        if os.path.exists(os.path.join(path, "prefiberswap")):
            print(path, "already run")
            continue
        plPlug = glob.glob(path + "/plPlug*.par")
        if not len(plPlug) == 1:
            print("no unique plPlugmatch", plPlug, "skipping")
            continue
        #scanPars = parsePlPlug(plPlug[0])
        try:
            scanPars = parsePlPlug(plPlug[0])
        except:
            print("failed to parse", plPlug[0])
            continue
        print("moving old scan data to saved directory")
        savedPath = os.path.join(path, "prefiberswap")
        os.makedirs(savedPath)
        tocopy = glob.glob(path + "/*.png") + glob.glob(path + "/*.par") + glob.glob(path+"/*.pkl")
        for copyme in tocopy:
            base, filename = os.path.split(copyme)
            newfile = os.path.join(savedPath, filename)
            shutil.copyfile(copyme, newfile)
        os.remove(plPlug[0])
        os.remove(os.path.join(path, "detectionList.pkl"))
        try:
            runMapper.resolve(
                scanDir = path,
                plateID = scanPars["plateId"],
                cartID = scanPars["cartridgeId"],
                fscanID = scanPars["fscanId"],
                fscanMJD = scanPars["fscanMJD"],
                )
        except Exception as e:
            print(path, "FAILED!!!!")
            print(e)
            subprocess.call("cp %s/* %s"%(savedPath, path), shell=True)
            subprocess.call("rm -r %s"%savedPath, shell=True) 
            continue
        mjddir = "/data/mapper/%i"%scanPars["fscanMJD"]
        if not os.path.exists(mjddir):
            os.makedirs(mjddir)
        # copy the recently created plPlugMap to the mjddir
        plPlug = glob.glob(path+"/plPlug*.par")[0]
        print("copying", plPlug, "to", mjddir)
        subprocess.call("cp %s %s"%(plPlug, mjddir), shell=True)




        

