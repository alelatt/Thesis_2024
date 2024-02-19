import rebound
import time
import numpy as np
import matplotlib.pyplot as plt

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

def Extract_Ellipticity(sim):
	retarr = np.zeros(sim.N-1)

	for i in range(len(retarr)):
		retarr[i] = sim.particles[i+1].e

	return retarr


sim = rebound.Simulation()
sim.start_server(port=1234)
sim.units = ['mearth', 'day', 'AU']

sim.add('Sun')
sim.add(m = 1, a = 1, e = 0.7, Omega = 0)
sim.add(m = 10, a = 1.5, e = 0.8, Omega = 0.1)
sim.add(m = 100, a = 2, e = 0.95, Omega = np.pi/3)
sim.add(m = 100, a = 0.5, e = 0, Omega = 0)

sim.integrator = "mercurius"
sim.ri_mercurius.r_crit_hill = 5.
sim.ri_mercurius.L = "infinity"
sim.move_to_com()

ob1 = rebound.OrbitPlot(sim)
plt.show()
plt.close('all')

E0 = sim.energy()

times = np.linspace(0, 100000, 10000)
ecc_arr = np.zeros((sim.N-1, len(times)))
en_arr = np.zeros((sim.N, len(times)))

ecc_time = 0
en_time = 0

for i, timestep in enumerate(times):
	sim.integrate(timestep)

	t1 = time.time()
	ecc_arr[:,i] = Extract_Ellipticity(sim)
	t2 = time.time()
	en_arr[:,i] = Compute_Energy(sim)
	t3 = time.time()

	ecc_time += t2-t1
	en_time += t3-t2

E1 = sim.energy()

ob1 = rebound.OrbitPlot(sim)
plt.show()
plt.close('all')

if np.any(ecc_arr > 1):
	print("ECC > 1!")
if np.any(en_arr > 0):
	print("EN > 0!")

plt.figure("Eccentricity")
for i in range(len(ecc_arr)):
	plt.plot(ecc_arr[i,:])

plt.figure("Energy")
for i in range(len(en_arr)):
	plt.plot(en_arr[i,:])

plt.show()

print(ecc_time/10000, en_time/10000)

print(abs(E1-E0)/E0)