from kepler import *

orbit = createOrbitByState(StateVector(Vector3(7000, -200, 0), Vector3(0, 11, 0)), 398600.44, 0)

print(orbit.getStateByEpoch(10))