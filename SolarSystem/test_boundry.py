import numpy as np
import matplotlib.pyplot as plt

import os
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


def Simulation():
	sim = rebound.Simulation()

	sim.units = ['mearth', 'year', 'AU']

	sim.integrator = "mercurius"
	sim.ri_mercurius.r_crit_hill = 3.
	sim.ri_mercurius.L = "infinity"
	sim.collision = "direct"
	sim.collision_resolve = "halt"

	sim.add(m = 333000, r = 4.65e-3)
	sim.add(m = 3.18e2, P = 1, r = 4.78e-4)
	sim.add(m = 3.18e2, P = 2.005, r = 4.78e-4)
	sim.add(m = 3.18e2, P = 2.99, r = 4.78e-4)

	sim.start_server(port=1234)

	times = np.linspace(0, 1e6, int(1e6) + 1)
	flag = 0

	try:
		for time in times:
			try:
				sim.exit_min_distance = 0.1
				sim.move_to_com()
				sim.integrate(time, exact_finish_time=0)
				if np.any(Compute_Energy(sim) > 0) and flag == 0:
					sim.boundary = "open"
					maxa = 0
					maxP = 0
					for p in range(1, sim.N):
						if sim.particles[p].P > maxP:
							maxP = sim.particles[p].P
							maxa = sim.particles[p].a
					sim.configure_box(100*maxa)
					flag = 1
					print("Opened at {:.1f}".format(time))
			except rebound.Encounter:
				print("Near Miss", sim.t, time)
			finally:
				sim.exit_min_distance = 0
				sim.move_to_com()
				sim.integrate(time, exact_finish_time=0)
				print(sim.t, time)

	except rebound.Collision:
		print("Collision at {:.1f}".format(sim.t))
		return
	
	return


if __name__ == '__main__':
	Simulation()