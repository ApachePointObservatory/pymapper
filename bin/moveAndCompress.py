import os
import shutil

fromDirBase = "/data/savedLCOMapScans/"
toDirBase = "/data/rawmapper/"
walkDirs = os.walk(fromDirBase)
endDirs = {}
for path, dirs, files in walkDirs:
    if "fscan" in path:
        # keep only path after fromDir
        endDirs[path] = files
for fromDir, fileList in endDirs.iteritems():
    tailDir = path.split(fromDirBase)[-1]
    toDir = os.path.join(toDirBase, tailDir)
    if not os.path.exists(toDir):
        print("making directory", toDir)
        # os.makedirs(toDir)
    for file in fileList:
        newFile = os.path.exists(os.path.join(toDir, file))
        if not os.path.exists(newFile):
            print("copying file")
            fromFile = os.path.join(fromDir, file)
            toFile = os.path.join(toDir, file)
            # shutil.copyfile(fromFile, toFile)