import rebound
import time
import numpy as np
import matplotlib.pyplot as plt

sim = rebound.Simulation()
sim.start_server(port=1234)
sim.units = ['mearth', 'day', 'AU']

sim.add('Sun')
sim.add(m = 100, a = 1, e = 0, f = -0.391, Omega = 0)
sim.add(m = 0.001, a = 20, e = 0.95, f = -0.60, Omega = 0)

sim.integrator = "mercurius"
sim.ri_mercurius.r_crit_hill = 3.
sim.ri_mercurius.L = "infinity"
sim.move_to_com()

sim.status()
time.sleep(5)

ob1 = rebound.OrbitPlot(sim, particles=[1,2])
plt.show()
plt.close('all')

P = sim.particles[1].P

ecc_array = np.zeros(sim.N-1)

while True:
	sim.integrate(sim.t + P/10)

	for i in range(0, sim.N-1):
		ecc_array[i] = sim.particles[i+1].e
	
	if any(ecc_array >= 1) or sim.t >= 100000:
		print(ecc_array, sim.t, P)
		break

ob1 = rebound.OrbitPlot(sim, particles=[1,2])
plt.show()