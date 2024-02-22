import rebound
import time
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
import os


if __name__ == '__main__':
	sim = rebound.Simulation("problem.bin")
	sim.start_server(port=1234)
	sim.integrator = "whfast"
	#sim.ri_mercurius.r_crit_hill = 10.
	#sim.ri_mercurius.L = "C4"

	ob1 = rebound.OrbitPlot(sim)
	plt.show()
	plt.close('all')

	sim.integrate(-750000)

	ob1 = rebound.OrbitPlot(sim)
	plt.show()
	plt.close('all')

	sim.stop_server(port=1234)