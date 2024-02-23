import rebound
import time
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
import os


if __name__ == '__main__':
	sa = rebound.Simulationarchive("problem.bin")
	stoptime = sa[-1].t
	print(stoptime)

	sim = sa[0]
	sim.integrator = "mercurius"
	sim.ri_mercurius.r_crit_hill = 3.
	sim.ri_mercurius.L = "infinity"
	sim.move_to_com()

	sim.start_server(port=1234)

	ob1 = rebound.OrbitPlot(sim)
	plt.show()
	plt.close('all')

	times = np.linspace(0, stoptime, 1000)
	vel1 = np.zeros(len(times))
	dist1 = np.zeros(len(times))
	vel2 = np.zeros(len(times))
	dist2 = np.zeros(len(times))

	for i, time in enumerate(times):
		sim.integrate(time)
		vel1[i] = np.sqrt(sim.particles[-2].vx**2 + sim.particles[-2].vy**2 + sim.particles[-2].vz**2)
		dist1[i] = np.sqrt(sim.particles[-2].x**2 + sim.particles[-2].y**2 + sim.particles[-2].z**2)
		vel2[i] = np.sqrt(sim.particles[-1].vx**2 + sim.particles[-1].vy**2 + sim.particles[-1].vz**2)
		dist2[i] = np.sqrt(sim.particles[-1].x**2 + sim.particles[-1].y**2 + sim.particles[-1].z**2)

	ob1 = rebound.OrbitPlot(sim)

	plt.figure("Velocity")
	plt.plot(vel1, label = 'Planet 4')
	plt.plot(vel2, label = 'Planet 5')
	plt.legend(loc = 'best')
	plt.grid(True)

	plt.figure("Distance")
	plt.plot(dist1, label = 'Planet 4')
	plt.plot(dist2, label = 'Planet 5')
	plt.legend(loc = 'best')
	plt.grid(True)

	plt.show()
	plt.close('all')

	sim.stop_server(port=1234)