#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../utils.h"
#include "../../brute.h"
#include "../../stf.h"

float objective_function(float *p, void **extras); 

int main() {

	int ndim = 1;
	int nt = 10001;
	float dt = 0.002;
	int n_comb;
	int n;
	float tarr;
	int *n_par;
	float **par, **comb;
	void **extras;
	float *min;
	float *sv1;
	float pl[1], ph[1], pd[1];
	obj_fun_ptr func = &objective_function;
	float *full;
	int i,j;
	float t0;
	float slip, psv;
	float *tmp;
	float *stf, *t;
	int tmp_n;
	int t0_ind, tarr_ind;
	float Tr;
	float strike, dip, rake;
	FILE *fout, *fin;
	
	/////////////////////////////////////////////////
	// PARAMETER ESTIMATION SET UP
	/////////////////////////////////////////////////

	// allocate necessary stuff
	sv1 = (float*) malloc( sizeof(float) * nt );
	n_par = (int*) malloc( sizeof(int) * ndim );
	par = (float**) malloc( sizeof(float*) * ndim );
	extras = malloc( sizeof(void*) * 7 );

	// read test slip-rate function
	fin = fopen( "../scratch/sv1_1300_400_le.bin", "rb" );
	fread( sv1, sizeof(float), nt, fin );
	fclose( fin );

	// compute t0 and tarr
	t0_ind = find_first( sv1, nt, ">", 0.001 );

	// tarr is the first time after t0 where sv1 < 0.001
	tarr_ind = t0_ind + find_first( &sv1[t0_ind+1], nt-t0_ind-1, "<", 0.001 );

	// write these out to file in production
	t0 = t0_ind * dt;
	tarr = tarr_ind * dt;
	Tr = tarr - t0;
	printf("t0_ind=%d tarr_ind=%d t0=%f tarr=%f Tr=%f\n", t0_ind, tarr_ind, t0, tarr, Tr);


	// compute kinematic slip
	tmp_n = nt - (t0_ind + ((nt-1) - tarr_ind));
	slip = trapz( &sv1[t0_ind], tmp_n, dt );

	// compute peak slip velocity
	psv = maximum(sv1, nt);
	printf("k_slip=%f d_slip=%f psv=%f\n", slip, trapz(sv1, nt, dt), psv);

	// compute strike, dip, and rake
	

	// set-up extras
	extras[0] = (void*) sv1;
	extras[1] = (void*) &nt;
	extras[2] = (void*) &dt;
	extras[3] = (void*) &slip;
	extras[4] = (void*) &t0_ind;
	extras[5] = (void*) &tarr_ind;
	extras[6] = (void*) &Tr;

	// define parameters for estimation
	pl[0] = dt;
	ph[0] = Tr/2;
	pd[0] = dt;

	// create possibly jagged array of parameters		
	for (i=0; i<ndim; i++) 	{
		par[i] = arange(pl[i], ph[i], pd[i], &n_par[i]);
	}

	// generate all possible combinations of parameters
	comb = cartesian_product( par, n_par, &n_comb, ndim );
	full = (float*) malloc( sizeof(float) * n_comb);

	printf("n_comb=%d ndim=%d\n", n_comb, ndim);

	/////////////////////////////////////////////////
	// PARAMETER ESTIMATION
	/////////////////////////////////////////////////

	// check optimization
	min = fmin_brute( func, comb, ndim, n_comb, extras, full );
	t = arange( 0, nt*dt, dt, &tmp_n );
	stf = tinti( t, nt, Tr, min[0], t0, slip );

	// write estimated cost function
	fout = fopen( "../scratch/test_fdmin_cost_fun.bin", "wb" );
	fwrite( full, sizeof(float), n_comb, fout );
	fclose( fout );

	fout = fopen( "../scratch/test_fdmin_bf_tinti.bin", "wb" );
	fwrite( stf, sizeof(float), nt, fout );
	fclose( fout );

	printf("fmin(Ts) = (%f)\n", min[0]);

	// end nicely.
	dealloc2d_f(par, ndim);
	free(n_par);
	free(full);
	free(t);
	free(stf);
	return 0;

}

float objective_function(float *p, void **extras) {
	float *sv1;
	float *stf;
	float *t;
	float t0;
	int nt;
	int tmp_nt;
	int Tr;
	float slip;
	float dt;
	float res = 0.0;
	int t0_ind, tarr_ind;
	int i;

	// unpackage extras
	sv1 = (float*) extras[0];
	nt = *(int*) extras[1];
	dt = *(float*) extras[2];
	slip = *(float*) extras[3];
	t0_ind = *(int*) extras[4];
	tarr_ind = *(int*) extras[5];
	Tr = *(float*) extras[6];

	// if (invo_c == 0) printf("nt=%d dt=%f slip=%f t0_ind=%d tarr_ind=%d\n", nt, dt, slip, t0_ind, tarr_ind);

	t0 = t0_ind * dt;

	t = arange(0, dt*nt, dt, &tmp_nt);
	if (tmp_nt != nt) {
		printf("ERROR: sim nt and stf nt not equivalent. aborting!\n");
		exit(-1);
	}

	// compute analytical tinti function with parameter
	stf = tinti(t, nt, Tr, p[0], t0, slip);

	// compute residual as sum of squares
	for (i=t0_ind; i<tarr_ind; i++) {
		res = res + powf((sv1[i] - stf[i]), 2);
	}
	free(t);
	free(stf);

	return res;
}