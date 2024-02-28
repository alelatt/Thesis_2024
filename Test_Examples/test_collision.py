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

if __name__ == '__main__':
	sa = rebound.Simulationarchive("./simarchive/archive12.bin")

	rebound.OrbitPlot(sa[-1])
	plt.show()
	plt.close('all')

	sim = sa[-2]
	sim.integrator = "mercurius"
	sim.ri_mercurius.r_crit_hill = 3.
	sim.ri_mercurius.L = "infinity"
	sim.start_server(port=1234)

	times = np.linspace(-523629, -523630, 1000)
	dist = np.zeros(1000)

	for i, time in enumerate(times):
		sim.integrate(time, exact_finish_time = 0)
		dist[i] = Distance(sim, 1, 2)

	ob1 = rebound.OrbitPlot(sim)

	plt.figure("Distance")
	plt.plot(dist)

	plt.show()
	plt.close('all')

	sim.stop_server(port=1234)