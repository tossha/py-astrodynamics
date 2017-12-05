from linalg import *

EARTH_MU = 398600.4415

class StateVector:

	def __init__(self, pos, vel):
		self.position = pos
		self.velocity = vel

class SimpleOrbit:
	def __init__(self, ecc, sma, aop, inc, mu, ta):
		self.mu = mu
		self.sma = sma
		self.ecc = ecc
		self.aop = aop
		self.inc = inc
		self.ta = ta
		self.elems = {
			'sma': sma,
			'rp':  sma * (1 - ecc),
			'ra':  sma * (1 + ecc),
			'ecc': ecc,
			'aop': aop,
			'inc': inc
		}

		self.p = self.sma * (1 - self.ecc ** 2)
		self.r = self.p / (1 + self.ecc * math.cos(self.ta))
		self.c = math.sqrt(self.mu * self.p)
		self.vm = self.c / self.r
		self.vr = math.sqrt(self.mu / self.p) * self.ecc * math.sin(self.ta)

		self.gradients = {
			'sma': {
				'vr': 2 * self.sma ** 2 * self.vr / self.mu,
				'vm': 2 * self.sma ** 2 * self.vm / self.mu,
				'vb': 0
			},
			'ecc': {
				'vr': self.c ** 2 * self.vr / self.mu ** 2 / self.ecc,
				'vm': self.c / self.mu ** 2 / self.ecc * (self.c * self.vm - self.mu * self.r / self.sma),
				'vb': 0
			},
			'inc': {
				'vr': 0,
				'vm': 0,
				'vb': math.cos(self.aop + self.ta) / self.vm
			},
		}
		self.gradients['rp'] = {
			'vr': self.gradients['sma']['vr'] * (1 - self.ecc) - self.sma * self.gradients['ecc']['vr'],
			'vm': self.gradients['sma']['vm'] * (1 - self.ecc) - self.sma * self.gradients['ecc']['vm'],
			'vb': 0
		}
		self.gradients['ra'] = {
			'vr': self.gradients['sma']['vr'] * (1 + self.ecc) + self.sma * self.gradients['ecc']['vr'],
			'vm': self.gradients['sma']['vm'] * (1 + self.ecc) + self.sma * self.gradients['ecc']['vm'],
			'vb': 0
		}


class Orbit:

	def __init__(self, ecc, sma, aop, inc, loan, anomaly, epoch, mu, isAnomalytrue = True):
		self.mu = mu
		self.sma = sma
		self.ecc = ecc
		self.inc = inc
		self.loan = loan
		self.aop = aop
		self.epoch = epoch

		self.isElliptic = self.ecc < 1

		if isAnomalytrue:
			self.m0 = self.getMeanAnomalyByTrueAnomaly(anomaly)
		else:
			self.m0 = anomaly

		self.updateMeanMotion()

	def setSma(self, sma):
		self.sma = sma
		return self
		
	def setEcc(self, ecc):
		self.ecc = ecc
		return self
		
	def setAop(self, aop):
		self.aop = aop
		return self
		
	def setInc(self, inc):
		self.inc = inc
		return self
		
	def setLoan(self, loan):
		self.loan = loan
		return self
		
	def setTa(self, ta):
		self.ta = ta
		return self
		
	def setEpoch(self, epoch):
		self.epoch = epoch
		return self

	def setMu(self, mu):
		self.mu = mu
		return self

	def addPrecession(self, r, j2, epoch):
		rate = -3/2 * r**2 * j2 * math.cos(self.inc) * self.meanMotion / self.sma**2 / (1 - self.ecc**2)**2
		self.m0 = self.getMeanAnomalyByEpoch(epoch)
		self.loan += rate * (epoch - self.epoch)
		self.epoch = epoch
		return self

	def updateMeanMotion(self):
		self.meanMotion = math.sqrt(self.mu / abs(self.sma)) / abs(self.sma);

	def getMeanAnomalyByEpoch(self, epoch):
		return self.m0 + self.meanMotion * (epoch - self.epoch)

	def getTrueAnomalyByEpoch(self, epoch):
		return self.getTrueAnomalyByMeanAnomaly(
			self.getMeanAnomalyByEpoch(epoch)
		)

	def getMeanAnomalyByEccentricAnomaly(self, ea):
		if self.isElliptic:
			return ea - self.ecc * math.sin(ea)
		else:
			return self.ecc * math.sinh(ea) - ea

	def getMeanAnomalyByTrueAnomaly(self, ta):
		return self.getMeanAnomalyByEccentricAnomaly(
			self.getEccentricAnomalyByTrueAnomaly(ta)
		)

	def getEccentricAnomalyByEpoch(self, epoch):
		return self.getEccentricAnomalyByMeanAnomaly(
			self.getMeanAnomalyByEpoch(epoch)
		)

	def getEccentricAnomalyByMeanAnomaly(self, ma):
		maxIter = 30
		delta = 0.00000001
		M = ma
		E = 0
		F = 0
		i = 0

		if self.isElliptic:
			M = M / (2 * math.pi)

			M = 2.0 * math.pi * (M - math.floor(M))

			E = M if (self.ecc < 0.8) else math.pi

			F = E - self.ecc * math.sin(M) - M

			while ((abs(F) > delta) and (i < maxIter)):
				E = E - F / (1.0 - self.ecc * math.cos(E))
				F = E - self.ecc * math.sin(E) - M
				i = i + 1
		else:
			E = (math.log(2 * (abs(M) + 1/3)) + 1) / self.ecc + (1 - 1 / self.ecc) * math.asinh(abs(M) / self.ecc)
			E *= math.copysign(1, M)

			F = self.ecc * math.sinh(E) - E - M

			while ((abs(F) > delta) and (i < maxIter)):
				E = E - F / (self.ecc * math.cosh(E) - 1)
				F = self.ecc * math.sinh(E) - E - M
				i = i + 1

		return E

	def getEccentricAnomalyByTrueAnomaly(self, ta):
		if self.isElliptic:
			cos = math.cos(ta)
			sin = math.sin(ta)
			cosE = (self.ecc + cos) / (1 + self.ecc * cos)
			sinE = math.sqrt(1 - self.ecc * self.ecc) * sin / (1 + self.ecc * cos)
			ang = math.acos(cosE)

			return ang if (sinE > 0) else (2 * math.pi - ang)
		else:
			return 2 * math.atanh(math.tan(ta / 2) / math.sqrt((self.ecc + 1) / (self.ecc - 1)))

	def getOwnCoordsByTrueAnomaly(ta):
		r = self.sma * (1 - self.ecc**2) / (1 + self.ecc * math.cos(ta))
		return Vector3(r * math.cos(ta), r * math.sin(ta), 0)

	# see http://microsat.sm.bmstu.ru/e-library/Ballistics/kepler.pdf
	def getStateByEccentricAnomaly(self, ea):
		cos = math.cos(ea)
		sin = math.sin(ea)
		koeff = math.sqrt(self.mu / self.sma) / (1 - self.ecc * cos)

		pos = Vector3(
			self.sma * (cos - self.ecc),
			self.sma * math.sqrt(1 - self.ecc * self.ecc) * sin,
			0
		)

		vel = Vector3(
			-koeff * sin,
			koeff * math.sqrt(1 - self.ecc * self.ecc) * cos,
			0
		)

		return StateVector(
			pos.rotateZ(self.aop).rotateX(self.inc).rotateZ(self.loan),
			vel.rotateZ(self.aop).rotateX(self.inc).rotateZ(self.loan)
		)

	def getStateByEpoch(self, epoch):
		if self.isElliptic:
			return self.getStateByEccentricAnomaly(self.getEccentricAnomalyByEpoch(epoch))

		ta = self.getTrueAnomalyByEpoch(epoch)
		pos = self.getOwnCoordsByTrueAnomaly(ta)
		flightPathAngle = math.atan(self.ecc * math.sin(ta) / (1 + self.ecc * math.cos(ta)))
		vel = Vector3(math.sqrt(self.mu * (2 / pos.mag - 1 / self.sma)), 0, 0).rotateZ(ta + math.pi / 2 - flightPathAngle)
		
		return StateVector(
			pos.rotateZ(self.aop).rotateX(self.inc).rotateZ(self.loan),
			vel.rotateZ(self.aop).rotateX(self.inc).rotateZ(self.loan)
		)


	def getState(self):
		return self.getStateByEccentricAnomaly(self.getEccentricAnomalyByTrueAnomaly(self.ta))

	def getOrbitalVelocityByEpoch(self, epoch):
		state = self.getStateByEpoch(epoch)
		orbitRadial = state.pos.unit()
		orbitNormal = state.position.cross(state.velocity).unit()
		vr = state.velocity.dot(orbitRadial) # radial
		vb = state.velocity.dot(orbitNormal) # normal
		vm = state.velocity.dot(orbitNormal.cross(orbitRadial)) #transversal

		return [vm, vr, vb]

def createOrbitByState(state, mu, epoch = 0):
	pos = state.position
	vel = state.velocity
	posMag = pos.mag

	angMomentum = pos.cross(vel)

	loan = math.atan2(angMomentum.x, -angMomentum.y)
	inc = math.atan2(math.sqrt(angMomentum.x**2 + angMomentum.y**2), angMomentum.z)
	sma = (mu * posMag) / (2 * mu - posMag * vel.mag**2)
	e = math.sqrt(1 - (angMomentum.mag**2 / (mu * sma)))

	p = pos.rotateZ(-loan).rotateX(-inc)
	u = math.atan2(p.y , p.x)

	radVel = pos.dot(vel)
	cos = (sma*(1 - e*e) / posMag - 1) / e
	ta = math.pi if cos < -1 else (0 if cos > 1 else math.acos(cos))
	
	if radVel < 0:
		ta = -ta

	aop = (u - ta) if ((u - ta) > 0) else 2 * math.pi + (u - ta);

	return Orbit(e, sma, aop, inc, loan, ta, epoch, mu)
