import numpy as np
import matplotlib.pyplot as plt

import os
import rebound
import time as tm

from multiprocessing import Pool
from sim_library import *

from ctypes import cdll, c_double
clibheartbeat = cdll.LoadLibrary('/home/alelatt/Thesis_2024/Shared_Libs/heartbeat/heartbeat.so')

################## CONSTANTS ##################
orbit_fraction = 720
small_timestep = 1/8760 #fraction of timestep to use while events are ongoing
exit_enc_const = 2 #number of mutual hill radii below which a close encounter is detected
exit_esc_const = 1e-3 #U(planet, star)/T(star) under which an escape is detected

AU = 149597870700
REarth = 6371000/AU
RJup = 69911000/AU
RSun = 695700000/AU
MSun = 1.047348644e3	#in jupiter masses

start_time = 0
end_time = -int(1e4)
step_time = int(1e4)
step_fraction = 1000
N_processes = 4#12
plot_lines = 2#3
plot_columns = 2#4

def Parallel_Sim(niter):
	print(niter)
	sim = rebound.Simulation()

	sim.units = ['mjupiter', 'year', 'AU']

	sim.integrator = "ias15"
	sim.collision = "direct"
	sim.collision_resolve = "halt"
	sim.track_energy_offset = 1

	sim.heartbeat = clibheartbeat.heartbeat
	
	hashes = ["AU_Mic", "AU_Mic_b", "AU_Mic_c"]
	params = np.loadtxt("/home/alelatt/Thesis_2024/Test_Examples/AU_Mic/AU_Mic.txt")

	np.random.seed(niter)
	sim.add(m = 0.51*MSun, r = 0.744*RSun, hash = hashes[0])
	sim.add(m = params[0,9], r = params[0,12]*RJup, P = params[0,0]/365.25, e = np.random.uniform(params[0,4],params[0,5]), omega = params[0,6], hash = hashes[1])
	sim.add(m = params[1,9], r = params[1,12]*RJup, P = params[1,0]/365.25, e = np.random.uniform(params[1,4],params[1,5]), omega = params[1,6], hash = hashes[2])
	#sim.add(m = params[2,11], r = params[2,14]*REarth, P = params[2,0]/365.25, e = params[2,6], omega = np.random.uniform(params[2,6],params[2,7]), hash = hashes[3])

	t1 = tm.time()
	out_dir = Simulation(simulation = sim, t_start = start_time, t_end = end_time, step = step_time, step_fraction = step_fraction, hashes = hashes, orbit_fraction = orbit_fraction, small_timestep = small_timestep, exit_enc_const = exit_enc_const, exit_esc_const = exit_esc_const, info = True)

	return [out_dir, tm.time()-t1]



if __name__ == '__main__':
	folderpath = ""

	if folderpath == "":
		loctime = datetime.now()
		folderpath = "./{}_{}_{}_{}_{}_{}".format(loctime.year, loctime.month, loctime.day, loctime.hour, loctime.minute, loctime.second)

	if not os.path.exists(folderpath):
		os.mkdir(folderpath)
	os.chdir(folderpath)

	directories = []

	if len(os.listdir('./')) == 0:
		init_time = tm.time()
		pool = Pool()
		results = np.array(pool.map(Parallel_Sim, range(N_processes)))

		print("Total runtime {:.0f} min".format((tm.time() - init_time)/60))

		sing_runt = []

		for res_set in results:
			directories.append(res_set[0])
			print(res_set[1])
			sing_runt.append(float(res_set[1]))

		print("Average runtime of a process {:.0f} min".format(np.mean(np.array(sing_runt))/60))

	else:
		directories = [ f.name for f in os.scandir('./') if f.is_dir() ]

	plot_titles = ["Semimajor axis", "Eccentricity", "Inclination"]
	set_titles = ['tsteps', 'semiaxis', 'eccent', 'incl']
	output_sets = []
	for directory in directories:
		output_sets.append(np.load("{}/outputs.npz".format(directory)))

	for i in range(1, 4):
		fig, axs = plt.subplots(plot_lines, plot_columns, figsize=(10*plot_columns/plot_lines, 10), layout = 'tight')
		fig.suptitle(plot_titles[i-1])
		line = 0
		column = 0
		for j in range(N_processes):
			times = output_sets[j]['tsteps']
			values = output_sets[j][set_titles[i]]
			axs[line, column].plot(times, values, '.-')

			column += 1
			if column == plot_columns:
				column = 0
				line += 1

	plt.show()
