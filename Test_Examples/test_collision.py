import rebound
import time
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':
	sa = rebound.Simulationarchive("./simarchive/archive12.bin")
	print(sa[-1].t, Compute_Energy(sa[-1]))

	rebound.OrbitPlot(sa[-1])
	plt.show()
	plt.close('all')

	sim = sa[-2]
	sim.integrator = "mercurius"
	sim.ri_mercurius.r_crit_hill = 3.
	sim.ri_mercurius.L = "infinity"
	sim.start_server(port=1234)

	times = np.linspace(-523629, -523630, 10)

	for time in times:
		sim.integrate(time, exact_finish_time = 0)
		ob1 = rebound.OrbitPlot(sim)
		plt.show()
		plt.close('all')

	ob1 = rebound.OrbitPlot(sim)
	plt.show()
	plt.close('all')

	sim.stop_server(port=1234)