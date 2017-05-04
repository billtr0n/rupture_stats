#include <stdlib.h>
#include <stdio.h>
#include "../../utils.h"
#include "../../brute.h"

float objective_function(float *p, float **extras, int *ne); 

int main() {
	obj_fun_ptr func = &objective_function;
	float *full;
	int i,j;
	int ndim = 2;
	int n_comb;
	int *n_par;
	float **par, **comb;
	float *min;
	float pl[2], ph[2], pd[2];
	FILE *fout;

	// starting parameters
	pl[0] = -5.12;
	pl[1] = -5.12;
	ph[0] = 5.125;
	ph[1] = 5.125;
	pd[0] = 0.01;
	pd[1] = 0.01;

	// allocate necessary stuff
	n_par = (int*) malloc( sizeof(int) * ndim );
	par = (float**) malloc( sizeof(float*) * ndim );

	// create possibly jagged array of parameters		
	for (i=0; i<ndim; i++) {
		par[i] = arange(pl[i], ph[i], pd[i], &n_par[i]);
	}
	// generate all possible combinations of parameters
	comb = cartesian_product( par, n_par, &n_comb, ndim );
	full = (float*) malloc( sizeof(float) * n_comb);

	printf("n_comb=%d ndim=%d\n", n_comb, ndim);

	// check optimization
	min = fmin_brute( func, comb, ndim, n_comb, NULL, NULL, full );

	// write estimated cost function
	fout = fopen( "../scratch/test_fdmin_output.bin", "wb" );
	fwrite( full, sizeof(float), n_comb, fout );
	fclose( fout );

	printf("fmin(x,y) = (%f, %f)\n", min[0], min[1]);

	// end nicely.
	dealloc2d_f(par, ndim);
	free(n_par);
	free(full);
	return 0;

}

float objective_function(float *p, void **extras) {
	float r;
	r = rastrigin(p[0], p[1]);
	// r = ackley(p[0], p[1]);
	return r;
}

float empirical_objective_function(float *p, void **extras) {
	return 0.0;
}