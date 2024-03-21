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
	sim = rebound.Simulation()
	sim.units = ['mearth', 'year', 'AU']
	sim.integrator = "mercurius"
	sim.ri_mercurius.L = "infinity"
	sim.collision = "direct"
	sim.collision_resolve = "halt"
	sim.track_energy_offset = 1

	sim.start_server(port = 1234)

	np.random.seed(42)
	sim.add(m = 333000, r = 4.65e-3)
	sim.add(m = 3.18e2, P = np.random.uniform(1 - 0.01, 1 + 0.01), r = 4.78e-4)
	sim.add(m = 3.18e2, P = np.random.uniform(2 - 0.02, 2 + 0.02), r = 4.78e-4)
	sim.add(m = 3.18e2, P = np.random.uniform(3 - 0.03, 3 + 0.03), r = 4.78e-4)
	sim.add(m = 3.18e2, a = -1, e = 1.1, r = 4.78e-4)

	sim.move_to_com()

	minP = 100
	for i in range(1, sim.N):
		if sim.particles[i].P < minP:
			minP = sim.particles[i].P

	sim.dt = minP/360
	sim.ri_mercurius.r_crit_hill = 3.

	E0 = sim.energy()
	flag = 0
	lastN = sim.N

	times = np.arange(0, int(1e1) + 5, 5)

	try:
		for time in times:
			try:
				sim.usleep = 5e3

				if flag == 1:
					max_semi = 0
					for i in range(1,sim.N):
						if sim.particles[i].a > max_semi:
							max_semi = sim.particles[i].a
					sim.exit_max_distance = 25*max_semi
					print(sim.exit_max_distance, Distance(sim, 0, 4))

				sim.integrate(time, exact_finish_time=0)

			except rebound.Escape:
				print("Removed at", sim.t)
				index_rm = 1000
				for i in range(1,sim.N):
					if Distance(sim, 0, i) > sim.exit_max_distance:
						index_rm = i
				sim.remove(index = index_rm)
				flag = 0
				sim.move_to_com()
				sim.integrate(time, exact_finish_time=0)

			finally:
				if flag == 0 and np.any(Compute_Energy(sim) > 0):
					flag = 1
					print("Opened at {:.1f}".format(time))
	except:
		print("ERRORE")
	print(sim.t, sim.particles[1].x, sim.particles[2].x, sim.particles[3].x)
	return 1

if __name__ == '__main__':
	Nsym = 1

	t1 = time.time()
	pool = Pool()
	results = np.array(pool.map(Simulation, range(Nsym)))