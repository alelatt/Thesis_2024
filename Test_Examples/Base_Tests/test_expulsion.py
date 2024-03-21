import rebound
import time
import numpy as np
import matplotlib.pyplot as plt

def Distance(sim, index1, index2):
	pos1 = np.array([sim.particles[index1].x, sim.particles[index1].y, sim.particles[index1].z])
	pos2 = np.array([sim.particles[index2].x, sim.particles[index2].y, sim.particles[index2].z])
	return np.linalg.norm(pos2 - pos1)

def Kinetic(particle):
	return particle.m*(particle.vx**2 + particle.vy**2 + particle.vz**2)/2

def Potential(sim, index):
	U = 0
	for i in range(len(sim.particles)):
		if i != index:
			U += sim.G*sim.particles[i].m*sim.particles[index].m/Distance(sim, index, i)
	return U

def Compute_Energy(sim):
	energies = np.zeros(len(sim.particles))

	for i in range(len(energies)):
		energies[i] = Kinetic(sim.particles[i]) - Potential(sim, i)

	return energies


sim = rebound.Simulation()
sim.start_server(port=1234)
sim.units = ['mearth', 'day', 'AU']

sim.add('Sun')
sim.add(m = 100, a = 1, e = 0, f = -0.391, Omega = 0)
sim.add(m = 0.001, a = 20, e = 0.95, f = -0.60, Omega = 0)

sim.integrator = "mercurius"
sim.ri_mercurius.r_crit_hill = 3.
sim.ri_mercurius.L = "infinity"
sim.move_to_com()

#sim.status()
#time.sleep(5)

ob1 = rebound.OrbitPlot(sim, particles=[1,2])
plt.show()
plt.close('all')

P = sim.particles[1].P
step = P/360

ecc_array = np.zeros(sim.N-1)
en_array = np.zeros(sim.N)

while True:
	sim.integrate(sim.t + step)

	for i in range(0, sim.N-1):
		ecc_array[i] = sim.particles[i+1].e

	en_array = Compute_Energy(sim)

	if any(en_array > 0) or any(ecc_array >= 1) or sim.t >= 100000:
		print(ecc_array, sim.t, step)
		print(Compute_Energy(sim))
		break

ob1 = rebound.OrbitPlot(sim, particles=[1,2])
plt.show()