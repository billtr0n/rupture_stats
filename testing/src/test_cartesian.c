#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "../../utils.h"

int main() {
	
	float **A, **B;
	int *np, *c;
	int i, j;
    bool test;

	int ndim = 1;
	int ncom;

	float plow[1] = {1};
	float phigh[1] = {4};
	float pdel[1] = {1};

	// contains all total paramers
	A = (float**) malloc( ndim * sizeof(float*) );


	// counter initialized to zero
	c = zeros_i( ndim );

	// number of parameters 
	np = (int*) malloc( ndim * sizeof(int) );

	// create vector of parameters
	for (i=0; i<ndim; i++) {
		A[i] = arange( plow[i], phigh[i], pdel[i], &np[i] );
	}

	printf("Parameters\n");
	printf("==========\n");
	for (i=0; i<ndim; i++) {
		printf("{ ");
		for (j=0; j<np[i]; j++) {
			if (j == np[i]-1) {
				printf("%f", A[i][j]);
			} else {
				printf("%f, ", A[i][j]);
			}
		}
		printf(" }\n");
	}

	B =cartesian_product( A, np, &ncom, ndim );


	return 0;
}
