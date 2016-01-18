"""
Import this to do standard library and executable path setup for python programs.
This is not a program, it is in the bin directory so that is will be
automatically on the path.
"""
import sys, os
rootDir = os.path.normpath(os.path.dirname(os.path.dirname(sys.argv[0])))

# python path
pycbioDir = os.path.join(rootDir, "extern/pycbio")
if not os.path.exists(pycbioDir):
    raise Exception("can't fine pycbio directory: {}".format(pycbioDir))
sys.path = [pycbioDir, os.path.join(rootDir, "lib/py")] + sys.path

# executable PATH
clusterBin = "/cluster/bin/x86_64"
os.putenv("PATH", "{}:{}:{}".format(os.path.join(rootDir, "bin"),
                                    clusterBin,
                                    os.getenv("PATH")))
