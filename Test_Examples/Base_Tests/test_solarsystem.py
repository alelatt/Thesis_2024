import rebound
import time
import matplotlib.pyplot as plt

sim = rebound.Simulation()
sim.start_server(port=1234)
sim.units = ['mearth', 'day', 'AU']

sim.add('Sun')
sim.add('Mercury')
sim.add('Venus')
sim.add('399')
sim.add('301', primary = sim.particles[3])
sim.add('Mars')
sim.add('Jupiter')
sim.add('Saturn')
sim.add('Uranus')
sim.add('Neptune')

#sim.dt = 1
#sim.usleep = 0.0001
sim.integrator = 'whfast'
sim.ri_whfast.coordinates = "democraticheliocentric"
sim.move_to_com()

ob1 = rebound.OrbitPlot(sim, particles=[1,2,3,5,6,7,8,9])
ob2 = rebound.OrbitPlot(sim, particles=[4], primary = 3, fig = ob1.fig, ax = ob1.ax)
#ob1.xlim = [-5,5]
#ob1.ylim = [-5,5]
#ob1.update()
plt.show()
plt.close('all')

#sim.status()
sim.integrate(-4000)		#PER TORNARE INDIETRO NEL TEMPO BASTA IMPOSTARE UN LIMITE NEGATIVO

ob1 = rebound.OrbitPlot(sim, particles=[1,2,3,5,6,7,8,9])
ob2 = rebound.OrbitPlot(sim, particles=[4], primary = 3, fig = ob1.fig, ax = ob1.ax)
#ob1.xlim = [-5,5]
#ob1.ylim = [-5,5]
#ob1.update()
plt.show()