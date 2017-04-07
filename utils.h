/* creates 1d array of values starting with "start" ending with "end" and a
   length of "length". don't forget to free buffer in main when done! */
float* linspace(float start, float end, int length);

/* creates 1d array of values starting with "start" ending with "end" and a
   spacing of dx. don't forget to free buffer in main when done! */
float* arange(float start, float end, float dx, int *n);

/* creates zero initialized array.
    don't forget to free buffer in main when done! */
float* zeros(int n);

/* computes sum of 1d array */
float sum(float *buf, int n);

