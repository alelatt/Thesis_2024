import rebound
import time
import matplotlib.pyplot as plt

sim = rebound.Simulation()
sim.start_server(port=1234)
sim.units = ['mearth', 'day', 'AU']

time.sleep(5)

sim.add('Sun')
sim.add('Venus')
sim.add(m = 1, a = 1)
sim.add(m = 0.0123, a = 0.00257, primary = sim.particles[2])
sim.add('Mars')

#sim.dt = 1
sim.usleep = 100000
#sim.integrator = 'whfast'
sim.move_to_com()

ob1 = rebound.OrbitPlot(sim, particles=[1,2,4])
ob2 = rebound.OrbitPlot(sim, particles=[3], primary = 2, fig = ob1.fig, ax = ob1.ax)
ob1.xlim = [-5,5]
ob1.ylim = [-5,5]
ob1.update()
plt.show()
plt.close('all')

sim.status()
time.sleep(5)
sim.integrate(400)

ob1 = rebound.OrbitPlot(sim, particles=[1,2,4])
ob2 = rebound.OrbitPlot(sim, particles=[3], primary = 2, fig = ob1.fig, ax = ob1.ax)
ob1.xlim = [-5,5]
ob1.ylim = [-5,5]
ob1.update()
plt.show()