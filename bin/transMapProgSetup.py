"""
Import this to do standard library and executable path setup for python programs.
This is not a program, it is in the bin directory so that is will be
automatically on the path.
"""
import sys, os
rootDir = os.path.normpath(os.path.dirname(os.path.dirname(sys.argv[0])))

def _addExtern(module, relDir):
    modDir = os.path.join(rootDir,  "extern", module, relDir)
    if not os.path.exists(modDir):
        raise Exception("can't fine {} directory: {}".format(moduleDir, modDir))
    sys.path.insert(0, modDir)
_addExtern("pycbio", "lib")
_addExtern("pipettor", "build/lib")
sys.path.insert(0, os.path.join(rootDir,  "lib/py"))

# executable PATH
os.putenv("PATH", "{}:{}:{}".format(os.path.join(rootDir, "bin"),
                                    "/cluster/bin/x86_64",
                                    os.getenv("PATH")))
