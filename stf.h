/* returns tinti source time function based on analytical derivation 
   Note: all allocations done in main program.  the length of stf_out, should
   be equal to the length of t. */
void tinti(float* t, float Tr, float Tp, float t0, float slip, float* stf_out);

/* returns exponential source time function at times t. */
// float* exponential(float* t, float vpeak, float t0, float slip);

/* returns draeger source time function */
// float* draeger(float* t, float xi, float tau, float n, float slip);