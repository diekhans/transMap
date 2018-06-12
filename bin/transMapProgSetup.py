"""
Import this to do standard library and executable path setup for python programs.
This is not a program, it is in the bin directory so that is will be
automatically on the path.

Also provides some common initialization functions
"""
import sys
import os
import logging
rootDir = os.path.dirname(os.path.abspath(os.path.dirname(sys.argv[0])))

def _addExtern(module, relDir):
    modDir = os.path.join(rootDir, "extern", module, relDir)
    if not os.path.exists(modDir):
        raise Exception("can't find {} directory: {}".format(module, modDir))
    sys.path.insert(0, modDir)
_addExtern("pycbio", "lib")
_addExtern("pipettor", "build/lib")
sys.path.insert(0, os.path.join(rootDir,  "lib/py"))

# executable PATH, make sure we can override kent commands
os.environ["PATH"] = "{}:{}:{}:{}".format(os.path.join(rootDir, "bin"),
                                          os.path.expanduser("~/kent/bin/x86_64"),
                                          "/cluster/bin/x86_64",
                                          os.environ["PATH"])
from pycbio.sys import loggingOps
from pycbio.sys import fileOps
import pipettor

fileOps.setTmpEnv()

def addCmdOptions(parser):
    loggingOps.addCmdOptions(parser)

def setupFromCmd(opts):
    loggingOps.setupFromCmd(opts, sys.argv[0])
    pipettor.setDefaultLogging(logging.getLogger(), logging.DEBUG)
