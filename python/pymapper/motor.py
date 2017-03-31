""" slit head motor control
"""
from __future__ import division, absolute_import
import os

from twisted.internet.protocol import Protocol, ClientFactory
from twisted.internet.endpoints import TCP4ClientEndpoint
from twisted.internet import reactor
import numpy
# from twisted.internet.defer import Deferred

#@todo, implement timeouts
class MotorConfig(object):
    def __init__(self):
        self.configFile = os.path.join(os.getenv("PYMAPPER_DIR"), "etc", "motorConfig.dat")
        self.hostname = None
        self.port = None
        self.startPos = None
        self.endPos = None
        self.speed = None
        self.slitPos = None
        self.direction = None
        self.loadMe() # load from file and set attrs
        self.checkMe()


    def loadMe(self):
        slitPos = {}
        slitPopulate = False
        with open(self.configFile, "r") as f:
            lines = f.readlines()
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            line = line.lower()
            if line.startswith("hostname"):
                self.hostname = str(self.getlineValue(line))
            elif line.startswith("port"):
                self.port = int(self.getlineValue(line))
            elif line.startswith("startpos"):
                self.startPos = float(self.getlineValue(line))
            elif line.startswith("endpos"):
                self.endPos = float(self.getlineValue(line))
            elif line.startswith("speed"):
                self.speed = float(self.getlineValue(line))
            elif line.startswith("slitpos"):
                # begin populating dict
                slitPopulate = True
            elif line.startswith("}"):
                # done populating dict
                slitPopulate = False
                self.slitPos = slitPos
            elif slitPopulate:
                line = line.strip(",")
                fiber, motorPos = line.split(":")
                fiber = int(fiber)
                motorPos = float(motorPos)
                slitPos[fiber] = motorPos
        self.direction = numpy.sign(self.endPos-self.startPos)

    def getlineValue(self, line):
        return line.split("=")[-1].strip()

    def checkMe(self):
        if None in [
            self.hostname,
            self.port,
            self.startPos,
            self.endPos,
            self.speed,
            self.direction
            ]:
            raise RuntimeError("Some Missing Motor configuration")
        # check that all 300 fibers are in slit pos
        if numpy.array_equal(self.slitPos.keys(), range(1,301)) == 300:
            raise RuntimeError("Missing motor positions in slit pos config")

    def posFromTime(self, timestamp):
        """Return motor position for a given time
        """
        return self.startPos + self.direction*self.speed*timestamp

MOTOR_CONFIG = MotorConfig()

class Command(object):
    def __init__(self, cmdStr, callFunc=None, timeout=0):
        self.cmdStr = cmdStr
        self.callFuncs = []
        if callFunc is not None:
            self.callFuncs.append(callFunc)
        self.isDone = False

    def setDone(self):
        if self.isDone:
            raise RuntimeError("cannot set command %s done, already done!"%self.cmdStr)
        print("setting %s done!"%self.cmdStr)
        self.isDone = True
        for func in self.callFuncs:
            func()

    def addCallback(self, callFunc):
        self.callFuncs.append(callFunc)

class MotorProtocol(Protocol):

    def __init__(self, motorControllerInstance):
        self.mci = motorControllerInstance

    def dataReceived(self, data):
        """Called each time a line of data is received from the ASCII controller
        """
        self.mci.dataReceived(data)

    # def sendCommand(self, cmdStr):
    #     """Sent ascii text to the ASCII controller
    #     """
    #     self.transport.write("%s\n" % cmdStr)

    def connectionMade(self):
        """Called when a connection is made
        """
        print("connection made")

class MotorClientFactory(ClientFactory):
    def __init__(self, motorControllerInstance):
        self.mci = motorControllerInstance

    def startedConnecting(self, connector):
        print("Started to connect to motor.")

    def buildProtocol(self, addr):
        print("Connected to motor.")
        return MotorProtocol(self.mci)

    # def clientConnectionLost(self, connector, reason):
    #     print("Lost connection to motor.  Reason:", reason)

    # def clientConnectionFailed(self, connector, reason):
    #     print("Connection failed!")
        #raise RuntimeError("Connection to motor failed. Reason:%s"%reason)

# class MotorStatus(object):
#     def __init__(self):
#         self.speed = None
#         self.currentPosition = None
#         self.targetPosition = None
#         self.isHomed = None
#         self.laserOn = None

class MotorController(object):
    def __init__(self, readyCallback=None):
        """readyCallback called when MotorController is ready to scan!
        """
        self.readyCallback = readyCallback
        # self.status = MotorStatus()
        self.mcf = MotorClientFactory(self)
        self.protocol = None # will be set after connection made
        self.currCmd = Command(cmdStr="dummy")
        self.currCmd.setDone()
        self.commandQueue = []
        self.isHomed = False

    def addReadyCallback(self, readyCallback):
        self.readyCallback = readyCallback

    def connect(self):
        """Returns a deferred
        """
        point = TCP4ClientEndpoint(reactor, MOTOR_CONFIG.hostname, MOTOR_CONFIG.port)
        connDeferred = point.connect(self.mcf)
        connDeferred.addCallback(self.gotProtocol)
        # and then prepare the controller to scan!
        connDeferred.addCallback(self.prepareToScan)
        # if the connection failed, let us know
        connDeferred.addErrback(self.connFailed)

    def disconnect(self):
        print("disconnecting from ASCII server")
        return self.protocol.transport.loseConnection()
        print("killing twisted event loop")
        reactor.stop()

    def connFailed(self, failure):
        print("conn failed errback")
        print(str(failure))
        reactor.stop()
        # raise RuntimeError("conn failed", str(failure))

    def gotProtocol(self, protocol):
        self.protocol = protocol

    def prepareToScan(self, foo):
        print("preparing for scan")
        # foo is ignored arg passed via callback framework
        # could send a stop first...
        self.getStatus(callFunc=self.checkHomeThenMove)

    def scan(self, callFunc=None):
        print("beginning scan")
        self.move(MOTOR_CONFIG.endPos)
        self.laserOff(callFunc=callFunc)
        # send motor back to start position
        self.resetAfterScan()

    def resetAfterScan(self):
        print("resetAfterScan")
        # foo is ignored arg passed via callback framework
        # could send a stop first...
        # self.getStatus(callFunc=self.checkHomeThenMove)
        self.move(MOTOR_CONFIG.startPos, callFunc=self.disconnect)
        # try killin twisted event loop now?
        # self.protocol.transport.loseConnection()
        # reactor.stop()


    def checkHomeThenMove(self):
        if not self.isHomed:
            print("Slit Head Axis is not homed.  Home it before proceeding!")
            raise RuntimeError("Slit Head Axis is not homed.  Home it before proceeding!")
            reactor.stop()
            # raise RuntimeError("Slit Head Axis is not homed.  Home it before proceeding!")
        else:
            print("Axis is Homed!!")
            # move motor in position for scan.
            self.setSpeed(MOTOR_CONFIG.speed)
            self.move(MOTOR_CONFIG.startPos)
            self.laserOn(callFunc=self.readyCallback)

    def getStatus(self, callFunc=None):
        print("getStatus")
        return self.queueCommand("status", callFunc=callFunc)

    def setSpeed(self, value, callFunc=None):
        print("set speed to %.2f"%float(value))
        return self.queueCommand("speed %.2f"%float(value), callFunc=callFunc)

    def move(self, value, callFunc=None):
        print("move to %.2f"%float(value))
        return self.queueCommand("move %.2f"%float(value), callFunc=callFunc)

    def laserOn(self, callFunc=None):
        print("laser on")
        return self.queueCommand("lonn", callFunc=callFunc)

    def laserOff(self, callFunc=None):
        print("laser off")
        return self.queueCommand("loff", callFunc=callFunc)

    def dataReceived(self, data):
        if self.currCmd is None:
            print("unsolicited dataReceived: %s"%str(data))
            return # don't do anything with unsolicited output...
        for dataline in data.split("\n"):
            dataline = dataline.strip().lower()
            if not dataline:
                # ignore blank strings...
                continue
            print("laser output:", dataline)
            # right now I only care if the axis is homed
            # don't care about managing any other status bits,
            # however add a parser here to keep track of things
            # eg if status needs to be checked frequently...
            # data_lowered = data.lower()
            if "homed" in dataline:
                if "not_homed" in dataline:
                    self.isHomed = False
                else:
                    self.isHomed = True
            if dataline.endswith("ok"):
                # running command is done
                self.currCmd.setDone()

    def sendCommand(self, command):
        if not self.currCmd.isDone:
            raise RuntimeError("cannot send %s, currently busy with %s"%(command.cmdStr, self.currCmd.cmdStr))
        self.currCmd = command
        print("sending: ", command.cmdStr)
        self.protocol.transport.write(command.cmdStr)

    def queueCommand(self, cmdStr, callFunc=None):
        print("queueCommand", cmdStr)
        command = Command(cmdStr, callFunc=callFunc)
        command.addCallback(self.runQueue)
        self.commandQueue.append(command)
        self.runQueue()

    def runQueue(self):
        if not self.currCmd.isDone:
            # do nothing, command already executing
            return
        if self.commandQueue:
            # at least one command waiting to execute
            self.sendCommand(self.commandQueue.pop(0))




if __name__ == "__main__":
    mc = None
    def cleanup():
        global mc
        print("Cleaning up")
        mc.resetAfterScan()
    def imready():
        global mc
        print("I'm READY!!!!")
        mc.scan(cleanup)
    mc = MotorController(imready)
    # reactor.callLater(mc.resetAfterScan)
    reactor.run()

"""
status example:

SLIT_HEAD_AXIS:
__MOVE_ACTUAL_POSITION 0.0
__TARGET_POSITION 12.0000000
__DRIVE_STATUS: OFF
__MOTOR_CURRENT: 0.0
__DRIVE_SPEED_SP 0.89999998
__DRIVE_SPEED 0.89999998
__DRIVE_ACCEL 20
__DRIVE_DECEL 20
__MOVE_RANGE 0.0 - 155.000000
__HARDWARE_FAULT 0
__INSTRUCTION_FAULT 0
__HOMED
VERTICAL_AXIS:
__MOVE_ACTUAL_POSITION 13.1517000
__TARGET_POSITION 13.1999998
__DRIVE_STATUS: OFF
__MOTOR_CURRENT: 0.0
__DRIVE_SPEED_SP 50.0000000
__DRIVE_SPEED 50.0000000
__DRIVE_ACCEL 20
__DRIVE_DECEL 20
__MOVE_RANGE 0.0 - 950.000000
__HARDWARE_FAULT 0
__INSTRUCTION_FAULT 0
FOOT_SWITCH: OFF
LASER: OFF


not homed:

status
STATUS

SLIT_HEAD_AXIS:
__MOVE_ACTUAL_POSITION -0.01890000
__TARGET_POSITION 12.0000000
__DRIVE_STATUS: OFF
__MOTOR_CURRENT: 0.0
__DRIVE_SPEED_SP 1.00000000
__DRIVE_SPEED 1.00000000
__DRIVE_ACCEL 20
__DRIVE_DECEL 20
__MOVE_RANGE 0.0 - 155.000000
__HARDWARE_FAULT 0
__INSTRUCTION_FAULT 0
__NOT_HOMED
VERTICAL_AXIS:
__MOVE_ACTUAL_POSITION 13.1517000
__TARGET_POSITION 13.1999998
__DRIVE_STATUS: OFF
__MOTOR_CURRENT: 0.0
__DRIVE_SPEED_SP 50.0000000
__DRIVE_SPEED 50.0000000
__DRIVE_ACCEL 20
__DRIVE_DECEL 20
__MOVE_RANGE 0.0 - 950.000000
__HARDWARE_FAULT 0
__INSTRUCTION_FAULT 0
FOOT_SWITCH: OFF
LASER: OFF

OK


home
HOME

__SPEED: 1.00000000
__HOME_ACTUAL_POSITION 9.99999975e-05
OK
move 10
MOVE 10

__SPEED: 1.00000000
__MOVE_ACTUAL_POSITION 1.15330005
__MOVE_ACTUAL_POSITION 2.35339999
__MOVE_ACTUAL_POSITION 3.55539989
__MOVE_ACTUAL_POSITION 4.75740004
__MOVE_ACTUAL_POSITION 5.95730019
__MOVE_ACTUAL_POSITION 7.15939999
__MOVE_ACTUAL_POSITION 8.35939980
__MOVE_ACTUAL_POSITION 9.56140041
__MOVE_ACTUAL_POSITION 10.0000000
OK
home
HOME

__SPEED: 1.00000000
__HOME_ACTUAL_POSITION 9.12380028
__HOME_ACTUAL_POSITION 8.22379971
__HOME_ACTUAL_POSITION 7.32229996
__HOME_ACTUAL_POSITION 6.42070007
__HOME_ACTUAL_POSITION 5.52069998
__HOME_ACTUAL_POSITION 4.61920023
__HOME_ACTUAL_POSITION 3.71919990
__HOME_ACTUAL_POSITION 2.81769991
__HOME_ACTUAL_POSITION 1.91770005
__HOME_ACTUAL_POSITION 1.01619995
__HOME_ACTUAL_POSITION 0.11480000
__HOME_ACTUAL_POSITION 0.0
OK


move then stop

MOVE 10

__SPEED: 1.00000000
__MOVE_ACTUAL_POSITION 1.15530002
stop__MOVE_ACTUAL_POSITION 2.35739994

STOP


OK


ERROR INVALID COMMAND



status while move

MOVE 10

__SPEED: 1.00000000
status
STATUS

ERROR BUSY MOVING
__MOVE_ACTUAL_POSITION 1.15540004
__MOVE_ACTUAL_POSITION 2.35549998
__MOVE_ACTUAL_POSITION 3.55749989
__MOVE_ACTUAL_POSITION 4.75950003
__MOVE_ACTUAL_POSITION 5.95959997
__MOVE_ACTUAL_POSITION 7.15950012
__MOVE_ACTUAL_POSITION 8.36159992
__MOVE_ACTUAL_POSITION 9.56350040
__MOVE_ACTUAL_POSITION 10.0000000
OK


status while home

home
HOME

__SPEED: 1.00000000
status\
__HOME_ACTUAL_POSITION 9.12250042
STATUS\

ERROR BUSY HOMING
__HOME_ACTUAL_POSITION 8.22239971
__HOME_ACTUAL_POSITION 7.32240009
__HOME_ACTUAL_POSITION 6.42100000
__HOME_ACTUAL_POSITION 5.52099991
__HOME_ACTUAL_POSITION 4.61940002
__HOME_ACTUAL_POSITION 3.71790004
__HOME_ACTUAL_POSITION 2.81640005
__HOME_ACTUAL_POSITION 1.91649997
__HOME_ACTUAL_POSITION 1.01489997
__HOME_ACTUAL_POSITION 0.11340000
__HOME_ACTUAL_POSITION 0.0
OK

"""
