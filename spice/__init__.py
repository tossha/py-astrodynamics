from spiceypy import *
import os

commonKernelsPath = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'kernels' + os.sep

def loadKernel(name):
	if os.path.isfile(name):
		furnsh(name)
	else:
		commonPath = commonKernelsPath + name
		if os.path.isfile(commonPath):
			furnsh(commonPath)
		else:
			print('Kernel not found: ' + name)

# Leap seconds
loadKernel('naif0012.tls')

# Earth size and shape
loadKernel('pck00010.tpc')


loadKernel('gm_de431.tpc')

# Earth orientation (high accuracy)
loadKernel('earth_latest_high_prec.bpc')

# Earth orientation (lower accuracy)
loadKernel('earth_070425_370426_predict.bpc')

def getState(body, relativeTo, et, referenceFrame = 'ECLIPJ2000'):
	return list(spkezr(str(body), et, 'ECLIPJ2000', 'NONE', str(relativeTo))[0])
	
def getPosition(body, relativeTo, et, referenceFrame = 'ECLIPJ2000'):
	return getState(body, relativeTo, et, referenceFrame)[0:3]

def getVelocity(body, relativeTo, et, referenceFrame = 'ECLIPJ2000'):
	return getState(body, relativeTo, et, referenceFrame)[3:6]
