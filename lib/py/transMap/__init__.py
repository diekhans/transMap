"""standard setup for TransMap jobs"""
import os,sys,socket
from pycbio.sys import Pipeline

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
    pl = Pipeline.Procline(cmds, stdin=stdin, stdout=stdout)
    sys.stderr.write(pl.getDesc() + "\n")
    pl.wait()

