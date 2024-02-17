import rebound
import time
import numpy as np
import matplotlib.pyplot as plt

sim = rebound.Simulation()
sim.start_server(port=1234)
sim.units = ['mearth', 'day', 'AU']

sim.add('Sun')
sim.add(m = 30, a = 0.2, Omega = np.pi/2)
sim.add(m = 30, a = 0.4, Omega = 0)
sim.add(m = 1, a = 0.5, e = 0.8, Omega = np.pi/3)
sim.add(m = 1, a = 0.5, e = 0.8, Omega = np.pi)

#sim.dt = 1
#sim.usleep = 0.001
sim.integrator = 'whfast'
sim.ri_whfast.coordinates = "democraticheliocentric"
#sim.ri_whfast.min_dt = 1
sim.move_to_com()

ob1 = rebound.OrbitPlot(sim, particles=[1,2,3,4])
plt.show()
plt.close('all')

sim.status()
sim.integrate(100000)

ob1 = rebound.OrbitPlot(sim, particles=[1,2,3,4])
plt.show()