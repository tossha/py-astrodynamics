from kepler import *

orbit = createOrbitByState(StateVector(Vector3(7000, 0, 0), Vector3(0, 10.67, 0)), 398600.44, 0)

print(orbit.ecc)
print(orbit.sma)