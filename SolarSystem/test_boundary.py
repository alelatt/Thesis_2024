import numpy as np
import matplotlib.pyplot as plt

import os
import time
import rebound

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

def Update_Elems(sim, init_N, hashes, semiaxes, eccentricities, inclinations, timesteps):
	semi = np.zeros(init_N - 1)
	ecc = np.zeros(init_N - 1)
	incl = np.zeros(init_N - 1)
	for i in range(0, init_N - 1):
		try:
			semi[i] = sim.particles[hashes[i+1]].a
			ecc[i] = sim.particles[hashes[i+1]].e
			incl[i] = sim.particles[hashes[i+1]].inc
		except:
			semi[i] = np.nan
			ecc[i] = np.nan
			incl[i] = np.nan			

	semiaxes.append(semi.tolist())
	eccentricities.append(ecc.tolist())
	inclinations.append(incl.tolist())
	timesteps.append(sim.t)
	return semiaxes, eccentricities, inclinations, timesteps

def Handle_Expulsion(sim, time, hashes):
	index_rm = 1000
	for i in range(1,sim.N):
		if Distance(sim, 0, i) > sim.exit_max_distance:
			index_rm = i
	#sim.remove(index = index_rm)
	sim.boundary = "open"
	sim.configure_box(sim.exit_max_distance/2)
	sim.step()
	print("Removed Particle {:s}".format(hashes[index_rm]), sim.t, time)
	sim.boundary = "none"
	return

def Update_MaxDist(sim):
	max_semi = 0
	for i in range(1,sim.N):
		if sim.particles[i].a > max_semi:
			max_semi = sim.particles[i].a
	sim.exit_max_distance = 10*max_semi
	return

def Handle_Exceptions(sim, time, flag, init_N, hashes, semis, eccs, incls, tsteps, red_rate = 1e-3):
	'''
	flag = 1 : Close Encounter
	flag = 2 : Escape
	flag = 3 : Close Encounter and Escape
	'''
	flag_internal = flag

	timestep = 0
	if sim.t < time:
		timestep = min([sim.t + 1/8760, time])
	else:
		timestep = max([sim.t - 1/8760, time])

	exit_distance_min = sim.exit_min_distance
	exit_distance_max = sim.exit_max_distance

	if flag_internal == 1:
		print("Used flag 1 at {:.7f} before getting to {:.7f}, sub of {:.7f}".format(sim.t, timestep, time))
		sim.exit_min_distance = 0

		try:
			Update_MaxDist(sim)
			sim.integrate(timestep, exact_finish_time=0)

		except rebound.Escape:
			flag_internal = 3

		finally:
			semis, eccs, incls, tsteps = Update_Elems(sim, init_N, hashes, semis, eccs, incls, tsteps)

		sim.exit_min_distance = exit_distance_min

	if flag_internal == 2:
		index_rm = 1000
		for i in range(1,sim.N):
			if Distance(sim, 0, i) > sim.exit_max_distance:
				index_rm = i

		print("Used flag 2 because of {:s} at {:.7f} before getting to {:.7f}, sub of {:.7f}. Mass {:.3f}".format(hashes[index_rm], sim.t, timestep, time, sim.particles[index_rm].m))

		sim.exit_max_distance = np.inf

		try:
			if sim.particles[index_rm].m - red_rate >= 0:
				sim.particles[index_rm].m -= red_rate
			else:
				sim.particles[index_rm].m = 0

			if sim.particles[index_rm].m == 0:
				sim.boundary = "open"
				sim.configure_box(exit_distance_max/2)

			sim.integrate(timestep, exact_finish_time=0)

		except rebound.Encounter:
			flag_internal = 3

		finally:
			semis, eccs, incls, tsteps = Update_Elems(sim, init_N, hashes, semis, eccs, incls, tsteps)

		sim.exit_max_distance = exit_distance_max
		sim.boundary = "none"
	
	if flag_internal == 3:
		index_rm = 1000
		for i in range(1,sim.N):
			if Distance(sim, 0, i) > sim.exit_max_distance:
				index_rm = i

		print("Used flag 3 because of {:s} at {:.7f} before getting to {:.7f}, sub of {:.7f}".format(hashes[index_rm], sim.t, timestep, time))

		sim.exit_max_distance = np.inf
		sim.exit_min_distance = 0

		if sim.particles[index_rm].m - red_rate >= 0:
			sim.particles[index_rm].m -= red_rate
		else:
			sim.particles[index_rm].m = 0

		if sim.particles[index_rm].m == 0:
			sim.boundary = "open"
			sim.configure_box(exit_distance_max/2)
			print("{:s} deleted at {:.7f} before getting to {:.7f}, sub of {:.7f}".format(hashes[index_rm], sim.t, timestep, time))

		sim.integrate(timestep, exact_finish_time=0)

		semis, eccs, incls, tsteps = Update_Elems(sim, init_N, hashes, semis, eccs, incls, tsteps)

		sim.exit_min_distance = exit_distance_min
		sim.exit_max_distance = exit_distance_max
		sim.boundary = "none"
	return


def Simulation():
	sim = rebound.Simulation()

	sim.units = ['mearth', 'year', 'AU']

	sim.integrator = "mercurius"
	sim.ri_mercurius.r_crit_hill = 3.	#default 3
	sim.ri_mercurius.L = "infinity"
	#sim.ri_whfast.corrector = 17	#default 0
	#sim.ri_whfast.safe_mode = 0	#default 1
	#sim.ri_whfast.kernel = "lazy"	#default
	#sim.ri_ias15.epsilon = 1e-10	#default 1e-9
	sim.collision = "direct"
	sim.collision_resolve = "halt"
	sim.track_energy_offset = 1

	hashes = ["Sun", "Planet1", "Planet2", "Planet3", "Comet"]

	sim.add(m = 333000, r = 4.65e-3, hash = hashes[0])
	sim.add(m = 3.18e2, P = 1, r = 4.78e-4, hash = hashes[1])
	sim.add(m = 3.18e2, P = 2.01, r = 4.78e-4, hash = hashes[2])
	sim.add(m = 3.18e2, P = 2.99, r = 4.78e-4, hash = hashes[3])
	sim.add(m = 0.1, a = -0.5, e = 1.1, hash = hashes[4])

	sim.move_to_com()

	sim.start_server(port=1234)

	maxrhill = 0
	for i in range(1,sim.N):
		if sim.particles[i].rhill > maxrhill:
			maxrhill = sim.particles[i].rhill
	sim.exit_min_distance = 2*maxrhill

	minP = 100
	for i in range(1, sim.N):
		if sim.particles[i].P < minP:
			minP = sim.particles[i].P
	sim.dt = minP/360

	init_N = sim.N
	prev_N = init_N

	tsteps = []
	semis = []
	eccs = []
	incls = []

	E0 = sim.energy()

	times = np.arange(0, int(1e6) + 1, 1)

	try:
		for time in times:
			semis, eccs, incls, tsteps = Update_Elems(sim, init_N, hashes, semis, eccs, incls, tsteps)
			while abs(sim.t) < abs(time):
				try:
					sim.integrate(time, exact_finish_time=0)

				except rebound.Encounter:
					Handle_Exceptions(sim, time, 1, init_N, hashes, semis, eccs, incls, tsteps)
					'''
					#print("Near Miss", sim.t, time)
					sim.exit_min_distance = 0
					timestep = 0
					if sim.t < time:
						timestep = sim.t + 1/8760
					else:
						timestep = sim.t - 1/8760

					while abs(sim.t) < abs(timestep):
						try:
							Update_MaxDist(sim)
							sim.integrate(timestep, exact_finish_time=0)

						except rebound.Escape:
							Handle_Expulsion(sim, timestep, hashes)

						finally:
							semis, eccs, incls, tsteps = Update_Elems(sim, init_N, hashes, semis, eccs, incls, tsteps)

					sim.exit_min_distance = 2*maxrhill
					'''
				except rebound.Escape:
					Handle_Exceptions(sim, time, 2, init_N, hashes, semis, eccs, incls, tsteps)
					#Handle_Expulsion(sim, time, hashes)

	except rebound.Collision:
		print("Collision at {:.1f}".format(sim.t))
		print("DE/DT_integr", abs((sim.energy() + sim.energy_offset - E0)/E0)/sim.t)
		print("Energy offset", sim.energy_offset)
		print(sim.t, sim.particles[1].x, sim.particles[2].x, sim.particles[3].x)
		return np.array(semis), np.array(eccs), np.array(incls), np.array(tsteps)
	
	print(sim.t, sim.particles[1].x, sim.particles[2].x, sim.particles[3].x)
	return np.array(semis), np.array(eccs), np.array(incls), np.array(tsteps)


if __name__ == '__main__':
	t1 = time.time()
	semis, eccs, incls, tsteps = Simulation()
	print("Execution time {:.0f}".format(time.time() - t1))

	plt.figure("Semiaxes", figsize = [10,10])
	for i in range(len(semis[0,:])):
		plt.plot(tsteps, semis[:,i], '.-')

	plt.figure("Eccentricities", figsize = [10,10])
	for i in range(len(eccs[0,:])):
		plt.plot(tsteps, eccs[:,i], '.-')

	plt.figure("Inclinations", figsize = [10,10])
	for i in range(len(incls[0,:])):
		plt.plot(tsteps, incls[:,i], '.-')

	plt.show()