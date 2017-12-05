from spiceypy import *
import os

def loadKernel(name):
	if os.path.isfile(name):
		furnsh(name)
	else:
		commonPath = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'kernels' + os.sep + name
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
