import math
import numpy as np

class Vector3:

	def __init__(self, x, y, z):
		self.v = np.array([x, y, z])

	def __str__(self):
		return self.v.__str__()

	@property
	def x(self):
		return self[0]

	@x.setter
	def x(self, val):
		self[0] = val

	@property
	def y(self):
		return self[1]

	@y.setter
	def y(self, val):
		self[1] = val

	@property
	def z(self):
		return self[2]

	@z.setter
	def z(self, val):
		self[2] = val

	def __getitem__(self, key):
		return self.v[key]

	def __setitem__(self, key, value):
		self.v[key] = value
		return self

	def mag(self):
		return math.sqrt((self.v ** 2).sum())

	def rotateX(self, radians):
		sin = math.sin(radians)
		cos = math.cos(radians)
		return Vector3(
			self[0],
			self[1] * cos - self[2] * sin,
			self[1] * sin + self[2] * cos
		)

	def rotateY(self, radians):
		sin = math.sin(radians)
		cos = math.cos(radians)
		return Vector3(
		self[0] * cos + self[2] * sin,
		self[1],
		-self[0] * sin + self[2] * cos
		)

	def rotateZ(self, radians):
		sin = math.sin(radians)
		cos = math.cos(radians)
		return Vector3(
			self[0] * cos - self[1] * sin,
			self[0] * sin + self[1] * cos,
			self[2]
		)

	def cross(self, v):
		npRes = np.cross(self.v, v.v)
		return Vector3(npRes[0], npRes[1], npRes[2])

	def sub(self, v):
		npRes = np.subtract(self.v, v.v)
		return Vector3(npRes[0], npRes[1], npRes[2])

	def add(self, v):
		npRes = np.add(self.v, v.v)
		return Vector3(npRes[0], npRes[1], npRes[2])

	def dot(self, v):
		return np.dot(self.v, v.v)

	def unit():
		npRes = self.v / self.mag()
		return Vector3(npRes[0], npRes[1], npRes[2])

