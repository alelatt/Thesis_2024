import rebound
import time
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
import os

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
	t_init = time.time()
	sim = rebound.Simulation()
	sim.units = ['mearth', 'year', 'AU']
	sim.integrator = "mercurius"
	sim.ri_mercurius.L = "infinity"
	sim.collision = "direct"
	sim.collision_resolve = "halt"
	sim.track_energy_offset = 1

	np.random.seed(niter)
	sim.add(m = 333000, r = 4.65e-3)
	sim.add(m = 3.18e2, P = np.random.uniform(1 - 0.01, 1 + 0.01), r = 4.78e-4)
	sim.add(m = 3.18e2, P = np.random.uniform(2 - 0.02, 2 + 0.02), r = 4.78e-4)
	sim.add(m = 3.18e2, P = np.random.uniform(3 - 0.03, 3 + 0.03), r = 4.78e-4)

	sim.move_to_com()

	maxa = 0
	for i in range(1,sim.N):
		if sim.particles[i].a > maxa:
			maxrhill = sim.particles[i].a

	minP = 100
	for i in range(1, sim.N):
		if sim.particles[i].P < minP:
			minP = sim.particles[i].P

	sim.dt = minP/360
	sim.ri_mercurius.r_crit_hill = 3.
	#sim.ri_whfast.corrector = 17	#default 0
	#sim.ri_whfast.corrector2 =  1	#default 0
	#sim.ri_whfast.safe_mode = 0	#default 1
	#sim.ri_whfast.kernel = "lazy"	#default "default"
	#sim.ri_ias15.epsilon = 1e-9	#default 1e-9

	E0 = sim.energy()

	times = np.arange(0, int(1e6) + 1, 1)

	try:
		for timestep in times:
			sim.move_to_com()
			sim.integrate(timestep, exact_finish_time=0)

			if np.any(Compute_Energy(sim) > 0):
				print("System Unbound", niter)
				return sim.t, abs((sim.energy() + sim.energy_offset - E0)/E0)/sim.t, sim.t/(time.time() - t_init)
	except rebound.Collision:
		return sim.t, abs((sim.energy() + sim.energy_offset - E0)/E0)/sim.t, sim.t/(time.time() - t_init)

	return sim.t, abs((sim.energy() + sim.energy_offset - E0)/E0)/sim.t, sim.t/(time.time() - t_init)


if __name__ == '__main__':
	Nsym = 10

	t1 = time.time()
	pool = Pool()
	results = np.array(pool.map(Simulation, range(Nsym)))

	print("Total runtime %.d" %(time.time() - t1))
	print(results)
	print("Average DE/DT_sym", np.mean(results[:,1]))
	print("Average DT_sym/DT", np.mean(results[:,2]))