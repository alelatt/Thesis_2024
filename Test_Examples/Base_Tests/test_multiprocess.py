import rebound
import time
from multiprocessing import Pool
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

def Simulation(niter):
	print(niter)
	t_init = time.time()
	sim = rebound.Simulation()
	sim.units = ['mearth', 'year', 'AU']
	sim.integrator = "mercurius"
	sim.ri_mercurius.r_crit_hill = 3.
	sim.ri_mercurius.L = "infinity"

	sim.add(m = 333000)
	
	np.random.seed(niter)
	sim.add(m = 1, a = abs(np.random.normal(1, 0.2)), e = abs(np.random.normal(0.1, 0.05)), Omega = np.random.uniform(0, 2*np.pi))
	sim.add(m = 1, a = abs(np.random.normal(2, 0.2)), e = abs(np.random.normal(0.1, 0.05)), Omega = np.random.uniform(0, 2*np.pi))
	sim.add(m = 1, a = abs(np.random.normal(3, 0.2)), e = abs(np.random.normal(0.1, 0.05)), Omega = np.random.uniform(0, 2*np.pi))
	sim.add(m = 1, a = abs(np.random.normal(4, 0.2)), e = abs(np.random.normal(0.1, 0.05)), Omega = np.random.uniform(0, 2*np.pi))
	sim.add(m = 1, a = abs(np.random.normal(5, 0.2)), e = abs(np.random.normal(0.1, 0.05)), Omega = np.random.uniform(0, 2*np.pi))

	sim.move_to_com()

	times = np.linspace(0, -1e3, 100)

	for timestep in times:
		sim.integrate(timestep, exact_finish_time=0)

		if np.any(Compute_Energy(sim) > 0):
			return [1, time.time() - t_init, sim.t]

	return [0, time.time() - t_init, sim.t]


if __name__ == '__main__':
	t1 = time.time()
	pool = Pool()
	results = np.array(pool.map(Simulation, range(10)))
	
	print(results)
	print("Average time spent %.d, average sim time %.d" %(np.mean(results[:,1]), np.mean(results[:,2])))
	print("Stable %.d, unstable %.d" %(results[:,0].tolist().count(0), results[:,0].tolist().count(1)))
	print("Total runtime %.d" %(time.time() - t1))