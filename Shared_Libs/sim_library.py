import numpy as np
import matplotlib.pyplot as plt
import time as tm

import os
import rebound

from datetime import datetime

from ctypes import cdll, c_double
clibheartbeat = cdll.LoadLibrary('/home/alelatt/Thesis_2024/Shared_Libs/heartbeat/heartbeat.so')

def Distance(simulation, index1, index2):
	"""
	Computes distance between two particles in the simulation

	Inputs:
		simulation : rebound.Simulation
			Simulation object
		index1 : int
			Index of first particle
		index2 : int
			Index of second particle

	Outputs:
		float - Distance between the two particles
	"""
	part_1 = simulation.particles[index1]
	part_2 = simulation.particles[index2]
	pos1 = np.array([part_1.x, part_1.y, part_1.z])
	pos2 = np.array([part_2.x, part_2.y, part_2.z])
	return np.linalg.norm(pos2 - pos1)



def Compute_Distances(simulation):
	"""
	Computes distances between all particles

	Inputs:
		simulation : rebound.Simulation
			Simulation object
	
	Outputs:
		ndarray - Distances between all particles
	"""

	dist_array = []
	for i in range(1, simulation.N):
		j = 0
		while j<i:
			dist_array.append(Distance(simulation = simulation, index1 = i, index2 = j))
			j += 1
	return np.array(dist_array)



def Compute_Hill(simulation):
	"""
	Computes mutual Hill radius (if between planets) or Roche radius (if between
		planet and star) to be used to characterize close encounters.

	Inputs:
		simulation : rebound.Simulation
			Simulation object

	Outputs:
		ndarray - Hill or Roche radius between all particles
	"""

	hill_array = []
	for i in range(1, simulation.N):
		part_i = simulation.particles[i]
		j = 0
		while j<i:
			if j == 0:
				if (part_i.m == 0) and (Kinetic(simulation = simulation, index = i) - Potential(simulation = simulation, index = i) >= 0):
					hill_array.append(0)
				else:
					r_roche = part_i.r * ((2*simulation.particles[0].m/part_i.m)**(1./3.))
					hill_array.append(r_roche)
			else:
				part_j = simulation.particles[j]
				r_hill = ((part_i.a + part_j.a)/2) * (((part_i.m + part_j.m)/(3*simulation.particles[0].m))**(1./3.))
				hill_array.append(r_hill)
			j += 1
	return np.array(hill_array)



def Kinetic(simulation, index):
	"""
	Returns kinetic energy of a particle

	Inputs:
		simulation : rebound.Simulation
			Simulation object
		index : int
			Index of the particle
	"""
	part = simulation.particles[index]
	return part.m*(part.vx**2 + part.vy**2 + part.vz**2)/2



def Potential(simulation, index):
	"""
	Returns gravitational potential energy of a particle

	Inputs:
		simulation : rebound.Simulation
			Simulation object
		index : int
			Index of the particle
	"""

	U = 0
	for i in range(len(simulation.particles)):
		if i != index:
			U += simulation.G*simulation.particles[i].m*simulation.particles[index].m/Distance(simulation = simulation, index1 = index, index2 = i)
	return U



def Compute_Energy(simulation):
	"""
	Computes total energy of all particles

	Inputs:
		simulation : rebound.Simulation
			Simulation object

	Outputs:
		energies : ndarray
			Array of total energies
	"""

	energies = np.zeros(len(simulation.particles))
	for i in range(len(energies)):
		energies[i] = Kinetic(simulation = simulation, index = i) - Potential(simulation = simulation, index = i)
	return energies



def Update_Elems(simulation, init_N, hashes, outputs):
	"""
	Updates list of all outputs

	Inputs:
		simulation : rebound.Simulation
			Simulation object
		init_N : int
			Number of initial particles
		hashes : list
			List of strings with particle hashes
		outputs : list
			List of lists containing all outputs
	"""

	semi = np.zeros(init_N - 1)
	ecc = np.zeros(init_N - 1)
	incl = np.zeros(init_N - 1)

	for i in range(0, init_N - 1):
		try:
			orbit = simulation.particles[hashes[i+1]].orbit(primary = simulation.particles[0])
			semi[i] = orbit.a
			ecc[i] = orbit.e
			incl[i] = orbit.inc
		except:
			semi[i] = np.nan
			ecc[i] = np.nan
			incl[i] = np.nan

	outputs[0].append(simulation.t)
	outputs[1].append(semi.tolist())
	outputs[2].append(ecc.tolist())
	outputs[3].append(incl.tolist())
	return



def Update_Step(simulation, orbit_fraction):
	"""
	Updates starting guess for minimum timestep as fraction of shortest orbit

	Inputs:
		simulation : rebound.Simulation
			Simulation object
		orbit_fraction : float
			Fraction of shortest orbital period

	Outputs:
		float - Updated parameter	
	"""

	min_P = np.inf
	for i in range(1,simulation.N):
		if (simulation.particles[i].P < min_P) and (simulation.particles[i].P > 0):
			min_P = simulation.particles[i].P
	if simulation.dt > 0:
		return min_P/orbit_fraction
	else:
		return -min_P/orbit_fraction



def Output_to_file(simulation, outputs, fpath = ""):
	"""
	Saves simulation archive and zipped outputs

	Inputs:
		simulation : rebound.Simulation
			Simulation object
		outputs : list
			List of lists containing all outputs
		fpath : str
			Path for file saving
	"""

	if fpath == "":
		tm.sleep(np.random.uniform(0,1)/100)
		loctime = datetime.now()
		fpath = "./{}".format(loctime.microsecond)
		os.mkdir(fpath)
		return fpath

	if len(outputs[0]) != 0:
		simulation.save_to_file("{}/archive.bin".format(fpath))

		if os.path.isfile("{}/outputs.npz".format(fpath)):
			out_saved = np.load("{}/outputs.npz".format(fpath))
			np.savez_compressed("{}/outputs.npz".format(fpath),
				tsteps = np.concatenate((out_saved['tsteps'], np.array(outputs[0]))),
				semiaxis = np.concatenate((out_saved['semiaxis'], np.array(outputs[1]))),
				eccent = np.concatenate((out_saved['eccent'], np.array(outputs[2]))),
				incl = np.concatenate((out_saved['incl'], np.array(outputs[3]))))
		else:
			np.savez_compressed("{}/outputs.npz".format(fpath),
				tsteps = np.array(outputs[0]),
				semiaxis = np.array(outputs[1]),
				eccent = np.array(outputs[2]),
				incl = np.array(outputs[3]))

		if len(outputs[4]) > 0:
			fopen = open(fpath+"/register.txt", 'a+')
			fopen.writelines(outputs[4])
			fopen.close()

		for i in range(len(outputs)):
			outputs[i].clear()
	return

def Handle_Exceptions(simulation, register, init_N, hashes, outputs, en_array, orbit_fraction, small_timestep, red_rate = 1e-3, info = False):
	"""
	Handles close encounters and escapes, triggered by exceptions.
	During the handling the timesteps are reduced to small_timestep.
	If a close encounter is being handled this function keeps working until
		the close encounter is present.
	If an escape is being handled this function reduces the mass of the object
		by red_rate at every small_timestep untill the object's mass is zero,
		at which point the object is deleted.
	The register is first filled when the exception is raised and is thus not
		empty when this function is called. When an event ends it's removed from
		the register. If during the handling of one exception another event is
		detected, it's relevant informations are inserted into the register.
		The function exits when the register is empty (thus all events ended).

	Inputs:
		simulation : rebound.Simulation
			Simulation object
		register : list
			Register for exceptions
		init_N : int
			Number of initial particles
		hashes : list
			List of strings with particle hashes
		outputs : list
			List of lists containing all outputs
		en_array : ndarray
			Contains current reference energy ([0]) and past energy losses ([1])
		small_timestep : float
			Small timestep for integrations during exceptions handling
		red_rate : float
			Mass subtracted per small_timestep when handling escapes
		info : bool
			If True prints events to console
	"""

	old_dt = simulation.dt
	timestep = 0
	if old_dt > 0:
		timestep = small_timestep
	else:
		timestep = -small_timestep

	enc_const = c_double.in_dll(clibheartbeat,"exit_enc_const").value
	esc_const = c_double.in_dll(clibheartbeat,"exit_esc_const").value

	while (len(register[0]) != 0) or (len(register[1]) != 0):
		if old_dt > 0:
			simulation.dt = np.min([small_timestep, abs(old_dt)])
		else:
			simulation.dt = -np.min([small_timestep, abs(old_dt)])

		if len(register[0]) != 0:
			c_double.in_dll(clibheartbeat,"exit_enc_const").value = 0
			old_dt = Update_Step(simulation = simulation, orbit_fraction = orbit_fraction)

		if len(register[1]) != 0:
			c_double.in_dll(clibheartbeat,"exit_esc_const").value = 0
			for i in range(1,simulation.N):
				if simulation.G*simulation.particles[0].m*simulation.particles[i].m/Distance(simulation = simulation, index1 = 0, index2 = i) < esc_const*Kinetic(simulation = simulation, index = 0):
					if hashes[i] not in register[1]:
						register[1].append(hashes[i])
						outstr = "{:s} (Mass {:.3f}) escaped at {:.7f}\n".format(hashes[i], simulation.particles[i].m, simulation.t)
						outputs[4].append(outstr)
						if info == True:
							print(outstr)

			for hash_rm in register[1]:
				if simulation.particles[hash_rm].m == 0:
					simulation.remove(hash = hash_rm)
					register[1].remove(hash_rm)
					outstr = "{:s} deleted at {:.7f}\n".format(hash_rm, simulation.t)
					outputs[4].append(outstr)
					if info == True:
						print(outstr)
					if not register[1]:
						en_array[0] = simulation.energy()
						c_double.in_dll(clibheartbeat,"exit_esc_const").value = esc_const
				elif simulation.particles[hash_rm].m - red_rate < 0:
					simulation.particles[hash_rm].m = 0
				else:
					simulation.particles[hash_rm].m -= red_rate

		integ_time = simulation.t + timestep

		try:
			simulation.integrate(integ_time)

		except rebound.Encounter:
			register[0].append(simulation.t)
			outstr = "Close Encounter at {:.7f} before getting to {:.7f}\n".format(simulation.t, integ_time)
			outputs[4].append(outstr)
			if info == True:
				print(outstr)
			c_double.in_dll(clibheartbeat,"exit_enc_const").value = 0
			simulation.integrate(integ_time)

		except rebound.Escape:
			c_double.in_dll(clibheartbeat,"exit_esc_const").value = 0
			for i in range(1,simulation.N):
				if simulation.G*simulation.particles[0].m*simulation.particles[i].m/Distance(simulation = simulation, index1 = 0, index2 = i) < esc_const*Kinetic(simulation = simulation, index = 0):
					if hashes[i] not in register[1]:
						register[1].append(hashes[i])
						outstr = "{:s} (Mass {:.3f}) escaped at {:.7f}\n".format(hashes[i], simulation.particles[i].m, simulation.t)
						outputs[4].append(outstr)
						if info == True:
							print(outstr)

		finally:
			Update_Elems(simulation = simulation, init_N = init_N, hashes = hashes, outputs = outputs)
			
		if (len(register[0]) != 0) and (np.all(Compute_Distances(simulation = simulation) > enc_const*Compute_Hill(simulation = simulation))):
			register[0].pop()
			c_double.in_dll(clibheartbeat,"exit_enc_const").value = enc_const

	simulation.dt = old_dt
	return

def Simulation(simulation, t_start, t_end, step, step_fraction, hashes, orbit_fraction, small_timestep, exit_enc_const, exit_esc_const, info = False):
	"""
	Runs the complete simulation

	Inputs:
		simulation : rebound.Simulation
			Simulation object
		t_start : float
			Start time of the simulation
		t_end : float
			End time of the simulation
		step : float
			Timestep for saving simulation archives
		step_fraction : float
			Fraction of a step to save orbital elements
		orbit_fraction : float
			Fraction of shortest orbital period
		small_timestep : float
			Small timestep for integrations during exceptions handling
		exit_enc_const : float
			Multiple of Hill or Roche radius to trigger close encounters
		exit_esc_const : float
			Value of (planet-star potential)/(star kinetic energy) under which
				an escape is triggered
		info : bool
			If True prints events to console

	Outputs:
		out_dir : str
			Directory in which data is saved

	
	The output list is divided in sublists such that:
		output[0] is a 1D list containing the reference times for the elements
		output[1] is a 2D list containing the semimajor axes at each time
		output[2] is a 2D list containing the eccentricities at each time
		output[3] is a 2D list containing the inclinations at each time
		output[4] is a 1D list containing a log of events

	The register is used in handling exceptions:
		register[0] is a list of all ongoing close encounters events
		register[1] is a list of all ongoing escape events
		The register is filled once an exception is called and is left empty by
			the Handle_Exceptions() function once all concurrent events end

	The Energy_array is used in keeping track of energy losses.
		Energy_array[0] is a float of the current reference energy
		Energy_array[1] is a float of the current energy loss
		Energy loss is computed at time t as (E(t0) - E(t))/(E(t0) * (t0-t))
			so that is a relative energy loss per unit time
		During an escape, since an object's mass is reduced to zero, the energy
			loss isn't counted
	"""
	simulation.move_to_com()

	time_dir = np.inf
	if t_start < t_end:
		time_dir = 1
	else:
		time_dir = -1

	c_double.in_dll(clibheartbeat,"exit_esc_const").value = exit_esc_const
	c_double.in_dll(clibheartbeat,"exit_enc_const").value = exit_enc_const

	simulation.dt = time_dir*simulation.dt
	simulation.dt = Update_Step(simulation, orbit_fraction)

	init_N = simulation.N

	output = [[], [], [], [], []]
	out_dir = Output_to_file(simulation = simulation, outputs = output)

	register = [[], []]

	E0 = simulation.energy()
	Delta_E = 0
	Energy_array = np.array([E0, Delta_E])

	times = np.arange(t_start + time_dir*step, t_end + time_dir*step, time_dir*step)
	Update_Elems(simulation = simulation, init_N = init_N, hashes = hashes, outputs = output)

	try:
		for time in times:
			Output_to_file(simulation = simulation, outputs = output, fpath = out_dir)

			substep = abs(time - simulation.t)/step_fraction
			subtimes = np.arange(simulation.t + time_dir*substep, time + time_dir*substep, time_dir*substep)
			for subtime in subtimes:
				while abs(simulation.t) < abs(subtime):
					try:
						simulation.integrate(subtime, exact_finish_time=0)
						Update_Elems(simulation = simulation, init_N = init_N, hashes = hashes, outputs = output)

					except rebound.Escape:
						Energy_array[1] = Energy_array[0] - simulation.energy()
						Energy_array[0] = simulation.energy()
						for i in range(1,simulation.N):
							if simulation.G*simulation.particles[0].m*simulation.particles[i].m/Distance(simulation = simulation, index1 = 0, index2 = i) < exit_esc_const*Kinetic(simulation = simulation, index = 0):
								if hashes[i] not in register[1]:
									register[1].append(hashes[i])
									outstr = "{:s} (Mass {:.3f}) escaped at {:.7f}\n".format(hashes[i], simulation.particles[i].m, simulation.t)
									output[4].append(outstr)
									if info == True:
										print(outstr)
						Handle_Exceptions(simulation = simulation, register = register, init_N = init_N, hashes = hashes, outputs = output, en_array = Energy_array, orbit_fraction = orbit_fraction, small_timestep = small_timestep, info = info)

					except rebound.Encounter:
						register[0].append(simulation.t)
						outstr = "Close Encounter at {:.7f} before getting to {:.7f}\n".format(simulation.t, subtime)
						output[4].append(outstr)
						if info == True:
							print(outstr)
						Handle_Exceptions(simulation = simulation, register = register, init_N = init_N, hashes = hashes, outputs = output, en_array = Energy_array, orbit_fraction = orbit_fraction, small_timestep = small_timestep, info = info)

	except rebound.Collision:
		outstr = "Collision at {:.7f}\nDE/(E0*DT_integr) = {}".format(simulation.t, (Energy_array[1] + Energy_array[0] - simulation.energy())/(simulation.t*E0))
		output[4].append(outstr)
		if info == True:
			print(outstr)
		Update_Elems(simulation = simulation, init_N = init_N, hashes = hashes, outputs = output)
		Output_to_file(simulation = simulation, outputs = output, fpath = out_dir)
		return out_dir
	
	outstr = "Integration over at {:.7f}\nDE/(E0*DT_integr) = {}".format(simulation.t, (Energy_array[1] + Energy_array[0] - simulation.energy())/(simulation.t*E0))
	output[4].append(outstr)
	if info == True:
		print(outstr)
	Update_Elems(simulation = simulation, init_N = init_N, hashes = hashes, outputs = output)
	Output_to_file(simulation = simulation, outputs = output, fpath = out_dir)
	return out_dir