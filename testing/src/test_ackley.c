#include <stdlib.h>
#include <stdio.h>
#include "../../utils.h"


int main() {

	int i,j;
	float *xi;
	float t;
	float **r;
	int n;
	FILE * out;

	out = fopen("../scratch/ackley.bin", "wb");

	xi = arange(-5.12,5.13,0.1,&n);

	t = ackley(0.0, 0.0);
	printf("f(0,0) = %f\n", t);

	r = (float**) malloc(sizeof(float*)*n);
	for (i=0; i<n; i++) {
		r[i] = (float*) malloc(sizeof(float)*n);
	}

	// printf( "n = %d\n", n );
	// for (i=0; i<n; i++) {
	// 	for (j=0; j<n; j++) {
	// 		r[i][j] = ackley(xi[i], xi[j]);
	// 	}
	// }

	// for (i=0; i<n; i++) {
	// 	fwrite(r[i], sizeof(float), n, out);
	// }

	free(xi);
	free(r);

	return 0;
}