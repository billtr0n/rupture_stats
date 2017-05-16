#include <stdio.h>
#include <stdlib.h>

#include "../../utils.h"
#include "../../stf.h"

int main() {
	float t0 = 1.0;
	float slip = 1.0;
	int i;
	float *t; 
	float *stf;
	int nt;
	float test;

	// test limiting conditions on function
	
	// define time vector
	t = arange(0, 10, 0.001, &nt);

	// generate tinti source-time function
	float Ts = 1.0;
	float Tr = 2.0;
	stf = tinti(t, nt, Tr, Ts, t0, slip);
	for (i=0; i<nt; i++) {
		printf("(%f) %f\n", t[i], stf[i]);
	}
    free(stf);

	 
	//   testing values at function boundaries.  
	//  note: values of inf are OK so long as every interval has an inequality 
	//   that provides finite values. 
	

	// test boundary conditions for all the integral ranges
	// printf("Case 1: Tr > 2*Ts\n");
	// printf("r1: %f\n", r1(0, Tr));
	// printf("r1: %f\n", r1(Ts, Tr));
	// printf("r2: %f\n", r2(Ts, Tr, Ts));
	// printf("r2: %f\n", r2(2 * Ts, Tr, Ts));
	// printf("r3: %f\n", r3(2 * Ts, Tr, Ts));
	// printf("r3: %f\n", r3(Tr, Tr, Ts));
	// printf("r4: %f\n", r4(Tr, Tr, Ts));
	// printf("r4: %f\n", r4(Tr + Ts, Tr, Ts));
	// printf("r5: %f\n", r5(Tr + Ts, Tr, Ts));
	// printf("r5: %f\n", r5(Tr + 2*Ts, Tr, Ts));

	// // second range only differs in the "third" integration range
	Ts = 1.0;
	Tr = 2.0;
	// generate tinti source-time function
	stf = tinti(t, nt, Tr, Ts, t0, slip);
	printf("Case 2: Ts < Tr < 2*Ts\n");
	printf("r1: %f\n", r1(0, Tr));
	printf("r1: %f\n", r1(Ts, Tr));
	printf("r2: %f\n", r2(Tr, Tr, Ts));
	printf("r2: %f\n", r2(Tr, Tr, Ts));
	printf("sr3: %f\n", sr3(Tr, Tr, Ts));
	printf("sr3: %f\n", sr3(2 * Ts, Tr, Ts));
	printf("r4: %f\n", r4(2 * Ts, Tr, Ts));
	printf("r4: %f\n", r4(Tr + Ts, Tr, Ts));
	printf("r5: %f\n", r5(Tr + Ts, Tr, Ts));
	printf("r5: %f\n", r5(Tr + 2*Ts, Tr, Ts));

	free(stf);
	free(t);

	return 0;
}
