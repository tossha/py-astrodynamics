from linalg import *

# EARTH_MU = 398600.4415
EARTH_MU = 398600.8

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

	def __init__(self, ecc, sma, aop, inc, epoch, mu, ta = False, loan = 0, ma = False):
		self.mu = mu
		self.sma = sma
		self.ecc = ecc
		self.inc = inc
		self.loan = loan
		self.aop = aop
		self.epoch = epoch
		self.ta = ta
		self.m0 = ma if ma != False else self.getMeanAnomalyByTrueAnomaly(ta)

		self.updateMeanMotion()

	# @property
	# def m0(self):
	# 	# if self.m0 != False:
	# 	# 	return self.m0
	# 	# self.m0 = self.getMeanAnomalyByTrueAnomaly(self.ta)
	# 	# return self.m0
	# 	return self.getMeanAnomalyByTrueAnomaly(self.ta)

	@property
	def p(self):
		return self.sma * (1 - self.ecc ** 2)

	@property
	def r(self):
		return self.p / (1 + self.ecc * math.cos(self.ta))

	@property
	def c(self):
		return math.sqrt(self.mu * self.p)

	@property
	def vm(self):
		return self.c / self.r

	@property
	def vr(self):
		return math.sqrt(self.mu / self.p) * self.ecc * math.sin(self.ta)

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

	def updateMeanMotion(self):
		self.meanMotion = math.sqrt(self.mu / self.sma) / self.sma;

	def getMeanAnomalyByEpoch(self, epoch):
		return self.m0 + self.meanMotion * (epoch - self.epoch)

	def getMeanAnomalyByEccentricAnomaly(self, ea):
		return ea - self.ecc * math.sin(ea)

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
		M = ma / (2.0 * math.pi)
		E = 0
		F = 0
		i = 0

		M = 2.0 * math.pi * (M - math.floor(M))

		E = M if (self.ecc < 0.8) else math.pi

		F = E - self.ecc * math.sin(M) - M

		while ((abs(F) > delta) and (i < maxIter)):
			E = E - F / (1.0 - self.ecc * math.cos(E))
			F = E - self.ecc * math.sin(E) - M
			i = i + 1

		return E

	def getEccentricAnomalyByTrueAnomaly(self, ta):
		cos = math.cos(ta)
		sin = math.sin(ta)
		cosE = (self.ecc + cos) / (1 + self.ecc * cos)
		sinE = math.sqrt(1 - self.ecc * self.ecc) * sin / (1 + self.ecc * cos)
		ang = math.acos(cosE)

		return ang if (sinE > 0) else (2 * math.pi - ang)

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
		return self.getStateByEccentricAnomaly(self.getEccentricAnomalyByEpoch(epoch))

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

# see http://microsat.sm.bmstu.ru/e-library/Ballistics/kepler.pdf
def createOrbitByState(state, mu, epoch = 0):
	pos = state.position;
	vel = state.velocity;
	posMag = pos.mag()

	angMomentum = pos.cross(vel);

	raan = math.atan2(angMomentum.x, -angMomentum.y);
	inc = math.atan2(math.sqrt(angMomentum.x**2 + angMomentum.y**2), angMomentum.z);
	sma = (mu * posMag) / (2.0 * mu - posMag * vel.mag()**2);
	e = math.sqrt(1.0 - (angMomentum.mag()**2 / (mu * sma)));

	p = pos.rotateZ(-raan).rotateX(-inc);
	u = math.atan2(p.y , p.x);

	radVel = pos.dot(vel) / posMag;
	cosE = (sma - posMag) / (sma * e);
	sinE = (posMag * radVel) / (e * math.sqrt(mu * sma));
	
	ta = math.atan2(math.sqrt(1.0 - e**2) * sinE, cosE - e);
	ta = ta if (ta > 0) else (ta + 2 * math.pi);

	aop = (u - ta) if ((u - ta) > 0) else 2 * math.pi + (u - ta);

	# print(e, sma, aop/math.pi*180, inc/math.pi*180)

	return Orbit(e, sma, aop, inc, epoch, mu, ta, raan)
