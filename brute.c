#include <stdlib.h>
#include <stdbool.h>
#include "utils.h"
#include "brute.h"
#include <stdio.h>

/* the object function will most likely be defined in another module or the main program.  this file contains only the actual optimization
   algorithm.  the empirical form also accepts a pointer to an array consisting of an empirical function to compute misfit against.
   see brute.h for parameter information.
*/
float* fmin_brute( obj_fun_ptr objective_function_ptr, float **p, int ndim, int np, void** extras, float* full_output ) {

	float **par;
	float *par_out;
	float **comb;
	float res, tmp;
	bool DEBUG = false;

	// indexing variables
	int i, pi, mi, di;

	// current minimum 
	float fmin = 0.0f;

	// check if NULL is passed and set flag accordingly
	bool full = (full_output != NULL) ? true : false;

	// allocate memory to store optimized parameters
	par_out = malloc(sizeof(float)*ndim);

	// main driver of the program
	for (i=0; i<np; i++) {

		// compute residual using cost function
		res = (*objective_function_ptr)(p[i], extras);

		// update or prime on first iteration
		if (i == 0 || res < fmin) {
			fmin = res;
			mi = i;
		}

		// save residual if full
		if (full) full_output[i] = res;
		if (DEBUG) printf("%f\n", res);
	}

	// deep copy of output
	for (i=0; i<ndim; i++) {
		par_out[i] = p[mi][i];
	}

	// return best fitting parameters
	return par_out;
}


/* use for unit testing 
int main() {

	return 0;
}*/
