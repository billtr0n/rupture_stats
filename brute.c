#include <stdlib.h>
#include <stdbool.h>
#include "utils.h"

/* the object function will most likely be defined in another module or the main program.  this file contains only the actual optimization
   algorithm.  the empirical form also accepts a pointer to an array consisting of an empirical function to compute misfit against.
   see brute.h for parameter information.
*/
float* fmin_brute( obj_fun_ptr func, int ndim, float *pl, float *ph, float *pd, float** extras, float* full_output ) {

	float **par;
	float **comb;
	int *n_par;
	float res, tmp;

	// indexing variables
	int i, pi, mi;

	// total number of combinations 
	int n_comb;

	// current minimum 
	float fmin = 0.0f;

	// check if NULL is passed and set flag accordingly
	bool full = (full_output != NULL) ? true : false;

	// create parameters
	n_par = (int*) malloc( sizeof(int) * ndim );
	par = (float**) malloc( sizeof(float*) * ndim );
	for (i=0; i<ndim; i++) {
		par[i] = arange(plow[i], phigh[i], pdel[i], &n_par[i]);
	}

	// generate all possible combinations of parameters
	comb = cartesian_product( par, n_par, &n_comb, ndim );

	// allocate memory to store global cost function
	if (full) {
		full_output = (float*) malloc( sizeof(float) * n_comb);
	}

	// main driver of the program
	for (i=0; i<n_comb; i++) {

		// compute residual using cost function
		res = func(comb[i], extras);

		// update or prime on first iteration
		if (i == 0 || res < fmin) {
			fmin = res;
			mi = i;
		}
 		
		// save residual if full
		if (full) full_output[i] = res;
	}

	// free stuff
	free(n_par);
	dealloc2d_f(par, ndim);
	dealloc2d_f(comb, n_comb);

	// return best fitting parameters
	return comb[mi];
}


/* use for unit testing 
int main() {

	return 0;
}*/