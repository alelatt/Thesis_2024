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
	print(niter)
	t_init = time.time()
	sim = rebound.Simulation()
	sim.units = ['mearth', 'year', 'AU']
	sim.integrator = "mercurius"
	sim.ri_mercurius.r_crit_hill = 3.
	sim.ri_mercurius.L = "infinity"

	sim.add(m = 333000)

	if niter == 0:
		sim.add(m = 3.18e2, P = 1)
		sim.add(m = 3.18e2, P = 2)
		sim.add(m = 3.18e2, P = 3)
	else:
		np.random.seed(niter)
		sim.add(m = 3.18e2, P = np.random.uniform(1 - 0.01, 1 + 0.01))
		sim.add(m = 3.18e2, P = np.random.uniform(2 - 0.02, 2 + 0.02))
		sim.add(m = 3.18e2, P = np.random.uniform(3 - 0.03, 3 + 0.03))

	sim.move_to_com()

	times = np.linspace(0, 1e7, int(1e7)+1)
	semiaxes = []

	for i, timestep in enumerate(times):
		sim.integrate(timestep, exact_finish_time=0)
		axes = np.zeros(sim.N)
		axes[0] = timestep
		for j in range(1, sim.N):
			axes[j] = sim.particles[j].a
		semiaxes.append(axes.tolist())
		if i%int(1e5) == 0:
			sim.save_to_file("./simarchive_res/archive{0}.bin".format(niter))
		if np.any(Compute_Energy(sim) > 0):
			sim.save_to_file("./simarchive_res/archive{0}.bin".format(niter))
			np.save("./simarchive_res/semi{0}.npy".format(niter), np.array(semiaxes))
			return 1
	np.save("./simarchive_res/semi{0}.npy".format(niter), np.array(semiaxes))
	return 0


if __name__ == '__main__':
	Nsym = 10

	if not os.path.exists("./simarchive_res"):
		os.makedirs("./simarchive_res")

	if not os.path.isfile("./simarchive_res/archive{0}.bin".format(Nsym-1)):
		t1 = time.time()
		pool = Pool()
		results = np.array(pool.map(Simulation, range(Nsym)))

		print(results)
		print("Total runtime %.d" %(time.time() - t1))
		print("Stable %.d, unstable %.d" %(results.tolist().count(0), results.tolist().count(1)))

	for i in range(Nsym):
		sa = rebound.Simulationarchive("./simarchive_res/archive{0}.bin".format(i))
		semiaxes = np.load("./simarchive_res/semi{0}.npy".format(i))
		plt.rcParams.update({'font.size': 18})
		title = "{:d}:P1 = {:.2f}, P2 = {:.2f}, P3 = {:.2f}".format(i,sa[0].particles[1].P, sa[0].particles[2].P, sa[0].particles[3].P)
		plt.figure(title, figsize = [10,10])
		for j in range(1, sa[0].N):
			plt.plot(semiaxes[:,0], semiaxes[:,j])
			plt.xticks(rotation = 25)
			plt.xlabel("Year since beginning")
			plt.ylabel("P [yrs]")
			plt.grid()

	plt.show()