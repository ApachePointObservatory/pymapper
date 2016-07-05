#!/usr/bin/env python
import os

def extractValue(logLine):
    """Get the float value following the ":" in a line
    """
    return float(logLine.split(":")[-1].strip())

def getScanParams(logfile):
    """Parse the logfile to determine the scan params
    """
    speed = None
    start = None
    end = None
    with open(logfile, "r") as f:
        logLines = f.readlines()
    for line in logLines:
        if "motor start pos" in line:
            start = extractValue(line)
        elif "motor end pos" in line:
            end = extractValue(line)
        elif "motor scan speed" in line:
            speed = extractValue(line)
        if not None in [start, end, speed]:
            break
    if None in [start, end, speed]:
        raise RuntimeError("Could not extract start, end, and/or speed: (%s, %s, %s)"
            %(str(start), str(end), str(speed)))
    return start, end, speed


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Determine fiber positions on the slit."
        )
    parser.add_argument("scanDir", help="""Directory containing scan""")
    scanDir = os.path.abspath(args.scanDir)
    if not os.path.exists(scanDir):
        raise RuntimeError("Scan directory does not exit: %s"%scanDir)
    logfile = os.path.join(scanDir, "scan.log")
    if not os.path.exists(logfile):
        raise RuntimeError("Could not locate scan log: %s"%logfile)
    start, end, speed = getScanParams(logfile)
    # get the ordered detections
    detectionListFile = os.path.join(scanDir, "detectionList.pkl")
    if not os.path.exists(detectionListFile):
        raise RuntimeError("Could not locate detection list file: %s"%(detectionListFile))
    pkl = open(detectionListFile, "rb")
    detectedFiberList = pickle.load(pkl)
    pkl.close()
    if not len(detectedFiberList)==300:
        raise RuntimeError("Number of detections %i != 300"%len(detectedFiberList))
    import pdb; pdb.set_trace()

if __name__ == "__main__":
    main()
