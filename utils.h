#include <stdbool.h>

/* creates 1d array of values starting with "start" ending with "end" and a
   length of "length". don't forget to free buffer in main when done! */
float* linspace(float start, float end, int length);

/* creates 1d array of values starting with "start" ending with "end" and a
   spacing of dx. don't forget to free buffer in main when done! */
float* arange(float start, float end, float dx, int *n);

/* creates zero initialized array of floats.
    don't forget to free buffer in main when done! */
float* zeros_f(int n);

/* creates zero initialized array of ints.
    don't forget to free buffer in main when done! */
int* zeros_i(int n);

/* computes sum of 1d array */
float sum(float *buf, int n);

/*
	computes cartesian product of jagged array A and returns all combinations in B. 

	inputs
		float **A : array containing all potential values of the parameters.
		int *np : length of A[i][]
		int ndim : number of dimensions
		int *ncom (by ref) : number of possible combinations


	returns
		float **B : size B[ncom][ndim] where ncom is the total number of combinations possible
*/
float** cartesian_product( float **A, int *np, int *ncom, int ndim );


/* 
	compares two floating point values to a value.  typically use eps = 1 / 2^20
*/
bool compare_f( float a, float b, float eps );

/* clean function to deallocate 2d float array */
void dealloc2d_f( float** buf, int n );

/* clean function to deallocate 2d int array */
void dealloc2d_i( int** buf, int n );