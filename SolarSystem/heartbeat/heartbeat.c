#include "rebound.h"
#include <math.h>
#include <stdlib.h>

/*
To compile:
move to heartbeat folder
gcc -c -O3 -fPIC heartbeat.c -o heartbeat.o
gcc -L. -shared heartbeat.o -o heartbeat.so -lrebound -Wl,-rpath='./heartbeat'
*/

/*
double max_distance = 0;

void heartbeat(struct reb_simulation* const r){
	const int N = r->N - r->N_var;
	for (int i=0;i<N;i++){
		double distance = reb_particle_distance(&r->particles[0], &r->particles[i]);
		if (distance>max_distance){
			double Kinetic = r->particles[i].m*(r->particles[i].vx*r->particles[i].vx + r->particles[i].vy*r->particles[i].vy + r->particles[i].vz*r->particles[i].vz)/2;
			double Potential = r->G*r->particles[0].m*r->particles[i].m/reb_particle_distance(&r->particles[0], &r->particles[i]);
			for (int j=1;j<N;j++){
				if (j != i){
					Potential += r->G*r->particles[i].m*r->particles[j].m/reb_particle_distance(&r->particles[i], &r->particles[j]);
				}
			}
			double Energy = Kinetic - Potential;
			if (Energy >= 0){
				r->status = REB_STATUS_ESCAPE;
			}
		}
	}
}
*/

double exit_esc_const = 0;
double exit_enc_const = 0;

void heartbeat(struct reb_simulation* const r){
	const int N = r->N - r->N_var;
	for (int i=0;i<N;i++){
		if (i > 0){
			double Kinetic = (r->particles[0].vx*r->particles[0].vx + r->particles[0].vy*r->particles[0].vy + r->particles[0].vz*r->particles[0].vz)/2;
			double Potential = r->G*r->particles[i].m/reb_particle_distance(&r->particles[0], &r->particles[i]);
			if (Potential < exit_esc_const*Kinetic){
				r->status = REB_STATUS_ESCAPE;
			}
		}
		for (int j=0;j<i;j++){
			if (j == 0){
				double r_roche = r->particles[i].r * pow(2*r->particles[0].m/r->particles[i].m, 1./3.);
				if (reb_particle_distance(&r->particles[0], &r->particles[i]) < exit_enc_const*r_roche){
					r->status = REB_STATUS_ENCOUNTER;
				}
			}
			else{
				struct reb_orbit o_i = reb_orbit_from_particle(r->G, r->particles[i], r->particles[0]);
				struct reb_orbit o_j = reb_orbit_from_particle(r->G, r->particles[j], r->particles[0]);
				double r_hill = ((abs(o_i.a) + abs(o_j.a))/2) * pow((r->particles[i].m + r->particles[j].m)/(3*r->particles[0].m), 1./3.);
				if (reb_particle_distance(&r->particles[i], &r->particles[j]) < exit_enc_const*r_hill){
					r->status = REB_STATUS_ENCOUNTER;
				}
			}
		}
	}
}