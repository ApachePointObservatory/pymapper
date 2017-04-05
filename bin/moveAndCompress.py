import os
import shutil
import subprocess
import glob

fromDirBase = "/home/lcomapper/scan/"
toDirBase = "/data/rawmapper/"
walkDirs = os.walk(fromDirBase)
endDirs = {}
for path, dirs, files in walkDirs:
    if "fscan" in path:
        tailDir = path.split(fromDirBase)[-1]
        toDir = os.path.join(toDirBase, tailDir)
        if not os.path.exists(toDir):
            print("creating", toDir)
            os.makedirs(toDir)
            for f in files:
                fromFile = os.path.join(path, f)
                toFile = os.path.join(toDir, f)
                if not os.path.exists(toFile):
                    print("creating file", toFile)
                    #shutil.copyfile(fromFile, toFile)
            # compress the files we just moved
            print("compressing fits files in %s"%toDir)
            #p = subprocess.Popen("fpack -D *.fits", cwd=toDir, shell=True)

#print("beginning to compress files via fpack")
#for path, dirs, files in os.walk(toDirBase):
#    if "fscan" in path:
#        if glob.glob(os.path.join(path, "*.fits")):
#            # uncompressed fits files exist
#            print("compressing images in ", path)
#            p = subprocess.Popen("fpack -D *.fits", cwd=path, shell=True)

