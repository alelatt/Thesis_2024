import rebound
import time
import matplotlib.pyplot as plt

sim = rebound.Simulation()
sim.start_server(port=1234)
sim.units = ['mearth', 'day', 'AU']

time.sleep(5)

sim.add('Sun')
sim.add('Earth')
sim.add('Moon', primary = 'Earth')

#sim.dt = 1
sim.usleep = 100000
#sim.integrator = 'whfast'
sim.move_to_com()

sim.status()
time.sleep(5)
sim.integrate(400)

rebound.OrbitPlot(sim)
plt.show()