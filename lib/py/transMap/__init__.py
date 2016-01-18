"""standard setup for TransMap jobs"""
import os,sys,socket
from pycbio.sys.pipeline import Procline

# FIXME: how many are still needed?

#def progSetup(

def getTmpExt():
    return "." + socket.gethostname() + "." + str(os.getpid()) + ".tmp"

def setTMPDIR():
    "setup TMPDIR env variable, also return dir path"
    if os.path.exists("/scratch/tmp"):
        tmpDir = "/scratch/tmp"
    else:
        tmpDir = "/var/tmp"
    os.environ["TMPDIR"] = tmpDir
    return tmpDir

def runCmds(cmds, stdin="/dev/null", stdout=None):
    """run a pipeline, printing command to stderr"""
    pl = Procline(cmds, stdin=stdin, stdout=stdout)
    sys.stderr.write(str(pl) + "\n")
    pl.wait()


def copyFile(src, dest):
    "copy a file with cp"
    runCmds(["cp", "-f", src, dest])

def getTwoBit(db):
    return "/hive/data/genomes/"+db +"/"+db+".2bit"

