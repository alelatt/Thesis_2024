import numpy as np
import matplotlib.pyplot as plt

import os
import rebound
import time as tm

from sim_library import *

from ctypes import cdll, c_double
clibheartbeat = cdll.LoadLibrary('/home/alelatt/Thesis_2024/Shared_Libs/heartbeat/heartbeat.so')

################## CONSTANTS ##################
orbit_fraction = 720
small_timestep = 1/8760 #one hour
exit_enc_const = 2 #number of mutual hill radii below which a close encounter is detected
exit_esc_const = 1e-3 #U(planet, star)/T(star) under which an escape is detected

if __name__ == '__main__':
	sim = rebound.Simulation()

	sim.units = ['mearth', 'year', 'AU']

	sim.integrator = "ias15"
	sim.collision = "direct"
	sim.collision_resolve = "halt"
	sim.track_energy_offset = 1

	sim.heartbeat = clibheartbeat.heartbeat
	
	hashes = ["Sun", "Planet1", "Planet2", "Planet3", "Comet"]
	sim.add(m = 333000, r = 4.65e-3, hash = hashes[0])
	sim.add(m = 3.18e2, P = 1, r = 4.78e-4, hash = hashes[1])
	sim.add(m = 3.18e2, P = 2.01, r = 4.78e-4, hash = hashes[2])
	sim.add(m = 3.18e2, P = 2.99, r = 4.78e-4, hash = hashes[3])
	sim.add(m = 3.18e1, a = -0.5, e = 1.1, hash = hashes[4])
	'''
	hashes = ["Sun", "Planet1", "Planet2", "Planet3"]
	sim.add(m = 333000, r = 4.65e-3, hash = hashes[0])
	sim.add(m = 0.81, P = 0.62, e = 0.007, r = 4.05e-5, hash = hashes[1])
	sim.add(m = 1, P = 1, e = 0.016, r = 4.26e-5, hash = hashes[2])
	sim.add(m = 0.1, P = 1.88, e = 0.09, r = 2.27e-5, hash = hashes[3])
	'''

	out_dir = []

	if len(out_dir) == 0:
		t1 = tm.time()
		out_dir = [Simulation(simulation = sim, t_start = 0, t_end = -int(1e4), step = int(1e3), hashes = hashes, orbit_fraction = orbit_fraction, small_timestep = small_timestep, exit_enc_const = exit_enc_const, exit_esc_const = exit_esc_const, info = True)]
		print("Execution time {:.0f}".format(tm.time() - t1))
	
	for directory in out_dir:
		outputs = np.load("{}/outputs.npz".format(directory))
		tsteps = outputs['tsteps']
		semis = outputs['semiaxis']
		eccs = outputs['eccent']
		incls = outputs['incl']

		plt.figure("Semiaxes {}".format(directory), figsize = [10,10])
		for i in range(len(semis[0,:])):
			plt.plot(tsteps, semis[:,i], '.-')

		plt.figure("Eccentricities {}".format(directory), figsize = [10,10])
		for i in range(len(eccs[0,:])):
			plt.plot(tsteps, eccs[:,i], '.-')

		plt.figure("Inclinations {}".format(directory), figsize = [10,10])
		for i in range(len(incls[0,:])):
			plt.plot(tsteps, incls[:,i], '.-')
	

	plt.show()