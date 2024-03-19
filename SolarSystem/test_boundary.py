import numpy as np
import matplotlib.pyplot as plt

import os
import time as tm
import rebound

from npy_append_array import NpyAppendArray
from ctypes import cdll, c_double
clibheartbeat = cdll.LoadLibrary('./heartbeat/heartbeat.so')

################## CONSTANTS ##################
orbit_fraction = 720
small_timestep = 1/8760 #one hour
exit_enc_const = 2 #number of mutual hill radii below which a close encounter is detected
exit_esc_const = 1e-3 #U(planet, star)/T(star) under which an escape is detected


def Distance(sim, index1, index2):
	pos1 = np.array([sim.particles[index1].x, sim.particles[index1].y, sim.particles[index1].z])
	pos2 = np.array([sim.particles[index2].x, sim.particles[index2].y, sim.particles[index2].z])
	return np.linalg.norm(pos2 - pos1)

def Compute_Distances(sim):
	dist_array = []
	for i in range(1, sim.N):
		j = 0
		while j<i:
			dist_array.append(Distance(sim, i, j))
			j += 1
	return np.array(dist_array)

def Kinetic(sim, index):
	return sim.particles[index].m*(sim.particles[index].vx**2 + sim.particles[index].vy**2 + sim.particles[index].vz**2)/2

def Potential(sim, index):
	U = 0
	for i in range(len(sim.particles)):
		if i != index:
			U += sim.G*sim.particles[i].m*sim.particles[index].m/Distance(sim, index, i)
	return U

def Compute_Energy(sim):
	energies = np.zeros(len(sim.particles))

	for i in range(len(energies)):
		energies[i] = Kinetic(sim, i) - Potential(sim, i)
	return energies

def Compute_Hill(sim):
	hill_array = []
	for i in range(1, sim.N):
		j = 0
		while j<i:
			if j == 0:
				if (sim.particles[i].m == 0) and (Kinetic(sim, i) - Potential(sim, i) >= 0):
					hill_array.append(0)
				else:
					r_roche = sim.particles[i].r * ((2*sim.particles[0].m/sim.particles[i].m)**(1./3.))
					hill_array.append(r_roche)
			else:
				r_hill = ((sim.particles[i].a + sim.particles[j].a)/2) * (((sim.particles[i].m + sim.particles[j].m)/(3*sim.particles[0].m))**(1./3.))
				hill_array.append(r_hill)
			j += 1
	return np.array(hill_array)

def Update_Elems(sim, init_N, hashes, outputs):
	semi = np.zeros(init_N - 1)
	ecc = np.zeros(init_N - 1)
	incl = np.zeros(init_N - 1)
	pot = 0
	for i in range(0, init_N - 1):
		try:
			semi[i] = sim.particles[hashes[i+1]].a
			ecc[i] = sim.particles[hashes[i+1]].e
			incl[i] = sim.particles[hashes[i+1]].inc
		except:
			semi[i] = np.nan
			ecc[i] = np.nan
			incl[i] = np.nan

	try:
		pot = sim.G*sim.particles[4].m*sim.particles[0].m/Distance(sim, 0, 4)
	except:
		pot = 0

	outputs[0].append(sim.t)
	outputs[1].append(semi.tolist())
	outputs[2].append(ecc.tolist())
	outputs[3].append(incl.tolist())
	outputs[5].append(Kinetic(sim, 0))
	outputs[6].append(pot)
	
	return

def Update_Params(sim):
	min_P = np.inf
	for i in range(1,sim.N):
		if (sim.particles[i].P < min_P) and (sim.particles[i].P > 0):
			min_P = sim.particles[i].P
	if sim.dt > 0:
		return min_P/orbit_fraction
	else:
		return -min_P/orbit_fraction

def Output_to_file(outputs, fpath = ""):
	if (not os.path.isdir(fpath)) and (fpath == ""):
		loctime = tm.localtime()
		fpath = "./{}_{}_{}_{}_{}".format(loctime.tm_year, loctime.tm_mon, loctime.tm_mday, loctime.tm_hour, loctime.tm_min,)
		os.mkdir(fpath)
		return fpath

	if len(outputs[0]) != 0:
		NpyAppendArray("{}/times.npy".format(fpath)).append(np.array(outputs[0]))
		NpyAppendArray("{}/elems.npy".format(fpath)).append(np.transpose(np.array(outputs[1:4]), (1, 0, 2)))
		NpyAppendArray("{}/kin.npy".format(fpath)).append(np.array(outputs[5]))
		NpyAppendArray("{}/pot.npy".format(fpath)).append(np.array(outputs[6]))
		file = open("{}/register.txt".format(fpath), 'a+')
		file.writelines(outputs[4])
		file.close()

		for i in range(len(outputs)):
			outputs[i].clear()

	return

def Handle_Exceptions(sim, register, init_N, hashes, outputs, en_array, red_rate = 1e-3):
	old_dt = sim.dt
	timestep = 0
	if old_dt > 0:
		timestep = small_timestep
	else:
		timestep = -small_timestep

	while (len(register[0]) != 0) or (len(register[1]) != 0):
		if old_dt > 0:
			sim.dt = np.min([small_timestep, abs(old_dt)])
		else:
			sim.dt = -np.min([small_timestep, abs(old_dt)])

		if len(register[0]) != 0:
			c_double.in_dll(clibheartbeat,"exit_enc_const").value = 0
			old_dt = Update_Params(sim)

		if len(register[1]) != 0:
			c_double.in_dll(clibheartbeat,"exit_esc_const").value = 0
			for i in range(1,sim.N):
				if sim.G*sim.particles[0].m*sim.particles[i].m/Distance(sim, 0, i) < exit_esc_const*Kinetic(sim, 0):
					if hashes[i] not in register[1]:
						register[1].append(hashes[i])
						outstr = "{:s} (Mass {:.3f}) escaped at {:.7f}\n".format(hashes[i], sim.particles[i].m, sim.t)
						outputs[4].append(outstr)
						print(outstr)

			for hash_rm in register[1]:
				if sim.particles[hash_rm].m == 0:
					sim.remove(hash = hash_rm)
					register[1].remove(hash_rm)
					outstr = "{:s} deleted at {:.7f}\n".format(hash_rm, sim.t)
					outputs[4].append(outstr)
					print(outstr)
					if not register[1]:
						en_array[0] = sim.energy()
						c_double.in_dll(clibheartbeat,"exit_esc_const").value = exit_esc_const
				elif sim.particles[hash_rm].m - red_rate < 0:
					sim.particles[hash_rm].m = 0
				else:
					sim.particles[hash_rm].m -= red_rate
				sim.ri_mercurius.recalculate_coordinates_this_timestep = 1

		integ_time = sim.t + timestep

		try:
			sim.integrate(integ_time)

		except rebound.Encounter:
			register[0].append(sim.t)
			outstr = "Close Encounter at {:.7f} before getting to {:.7f}\n".format(sim.t, integ_time)
			outputs[4].append(outstr)
			print(outstr)
			c_double.in_dll(clibheartbeat,"exit_enc_const").value = 0
			sim.integrate(integ_time)

		except rebound.Escape:
			c_double.in_dll(clibheartbeat,"exit_esc_const").value = 0
			for i in range(1,sim.N):
				if sim.G*sim.particles[0].m*sim.particles[i].m/Distance(sim, 0, i) < exit_esc_const*Kinetic(sim, 0):
					if hashes[i] not in register[1]:
						register[1].append(hashes[i])
						outstr = "{:s} (Mass {:.3f}) escaped at {:.7f}\n".format(hashes[i], sim.particles[i].m, sim.t)
						outputs[4].append(outstr)
						print(outstr)

		finally:
			Update_Elems(sim, init_N, hashes, outputs)
			
		if (len(register[0]) != 0) and (np.all(Compute_Distances(sim) > exit_enc_const*Compute_Hill(sim))):
			register[0].pop()
			c_double.in_dll(clibheartbeat,"exit_enc_const").value = exit_enc_const

	sim.dt = old_dt
	return

def Simulation(t_start, t_end, step):
	sim = rebound.Simulation()

	sim.units = ['mearth', 'year', 'AU']

	sim.integrator = "mercurius"
	#sim.ri_mercurius.r_crit_hill = 3.	#default 3
	#sim.ri_mercurius.L = "infinity"
	#sim.ri_mercurius.safe_mode = 0	#default 1
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

	sim.move_to_com()

	sim.start_server(port=1234)

	time_dir = np.inf
	if t_start < t_end:
		time_dir = 1
	else:
		time_dir = -1

	c_double.in_dll(clibheartbeat,"exit_esc_const").value = exit_esc_const
	c_double.in_dll(clibheartbeat,"exit_enc_const").value = exit_enc_const

	sim.dt = time_dir*sim.dt
	sim.dt = Update_Params(sim)

	rem_part = []
	init_N = sim.N

	tsteps = []
	semis = []
	eccs = []
	incls = []
	kin = []
	pot = []
	events = []
	output = [tsteps, semis, eccs, incls, events, kin, pot]
	out_dir = Output_to_file(output)

	closenc_register = []
	escape_register = []
	register = [closenc_register, escape_register]

	E0 = sim.energy()
	Delta_E = 0
	Energy_array = np.array([E0, Delta_E])

	times = np.arange(t_start + time_dir*step, t_end + time_dir*step, time_dir*step)
	Update_Elems(sim, init_N, hashes, output)

	try:
		for time in times:
			Output_to_file(output, out_dir)
			while abs(sim.t) < abs(time):
				try:
					sim.integrate(time, exact_finish_time=0)
					Update_Elems(sim, init_N, hashes, output)

				except rebound.Escape:
					Energy_array[1] = Energy_array[0] - sim.energy()
					Energy_array[0] = sim.energy()
					for i in range(1,sim.N):
						if sim.G*sim.particles[0].m*sim.particles[i].m/Distance(sim, 0, i) < exit_esc_const*Kinetic(sim, 0):
							if hashes[i] not in register[1]:
								register[1].append(hashes[i])
								outstr = "{:s} (Mass {:.3f}) escaped at {:.7f}\n".format(hashes[i], sim.particles[i].m, sim.t)
								output[4].append(outstr)
								print(outstr)
					Handle_Exceptions(sim, register, init_N, hashes, output, Energy_array)

				except rebound.Encounter:
					register[0].append(sim.t)
					outstr = "Close Encounter at {:.7f} before getting to {:.7f}\n".format(sim.t, time)
					output[4].append(outstr)
					print(outstr)
					Handle_Exceptions(sim, register, init_N, hashes, output, Energy_array)

	except rebound.Collision:
		outstr = "Collision at {:.7f}\nDE/(E0*DT_integr) = {}".format(sim.t, (Energy_array[1] + Energy_array[0] - sim.energy())/(sim.t*E0))
		output[4].append(outstr)
		print(outstr)
		Update_Elems(sim, init_N, hashes, output)
		Output_to_file(output, out_dir)
		print(sim.t, sim.particles[1].x, sim.particles[2].x, sim.particles[3].x)
		return out_dir
	
	outstr = "Integration over at {:.7f}\nDE/(E0*DT_integr) = {}".format(sim.t, (Energy_array[1] + Energy_array[0] - sim.energy())/(sim.t*E0))
	output[4].append(outstr)
	print(outstr)
	Update_Elems(sim, init_N, hashes, output)
	Output_to_file(output, out_dir)
	print(sim.t, sim.particles[1].x, sim.particles[2].x, sim.particles[3].x)
	return out_dir


if __name__ == '__main__':
	t1 = tm.time()
	out_dir = Simulation(0, int(1e5), 1000)
	print("Execution time {:.0f}".format(tm.time() - t1))
	'''
	plt.figure("Semiaxes 1", figsize = [10,10])
	for i in range(len(semis[0,:])):
		plt.plot(tsteps, semis[:,i], '.-')

	plt.figure("Eccentricities 1", figsize = [10,10])
	for i in range(len(eccs[0,:])):
		plt.plot(tsteps, eccs[:,i], '.-')

	plt.figure("Inclinations 1", figsize = [10,10])
	for i in range(len(incls[0,:])):
		plt.plot(tsteps, incls[:,i], '.-')

	plt.figure("Potential/Kinetic 1", figsize = [10,10])
	for i in range(len(incls[0,:])):
		plt.plot(tsteps, pot/kin, '.-')
	'''

	tsteps2 = np.load("{}/times.npy".format(out_dir))
	elems = np.transpose(np.load("{}/elems.npy".format(out_dir)), (1, 0, 2))
	kin2 = np.load("{}/kin.npy".format(out_dir))
	pot2 = np.load("{}/pot.npy".format(out_dir))
	semis2 = elems[0]
	eccs2 = elems[1]
	incls2 = elems[2]

	plt.figure("Semiaxes 2", figsize = [10,10])
	for i in range(len(semis2[0,:])):
		plt.plot(tsteps2, semis2[:,i], '.-')

	plt.figure("Eccentricities 2", figsize = [10,10])
	for i in range(len(eccs2[0,:])):
		plt.plot(tsteps2, eccs2[:,i], '.-')

	plt.figure("Inclinations 2", figsize = [10,10])
	for i in range(len(incls2[0,:])):
		plt.plot(tsteps2, incls2[:,i], '.-')

	plt.figure("Potential", figsize = [10,10])
	plt.plot(tsteps2, pot2, '.-')

	plt.figure("Kinetic", figsize = [10,10])
	plt.plot(tsteps2, kin2, '.-')

	plt.figure("Potential/Kinetic 2", figsize = [10,10])
	plt.plot(tsteps2, pot2/kin2, '.-')
	

	plt.show()