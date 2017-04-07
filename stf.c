#include <math.h>
#include "utils.h"

void* tinti(float* t, float Tr, float Ts, float t0, float slip, float* stf_out) {
	int i;
	int n;

	// get length of the time-series to compute the source time function
	n = length_f(t);

	// generate source time function
	for (i=0; i<n; i++) {
		// set of first ranges
		if (Tr > 2*Ts) {
			// 0 < t < ts
			if ((t >= 0) && (t < Ts)) {
				stf_out[i] = r1(t[i], Tr);
			}
			
			// ts < t < 2*ts
			else if ((t >= Ts) && (t < 2*Ts)) {
				stf_out[i] = r2(t[i], Tr, Ts);
			}

			// 2*ts < t < tr
			else if ((t>=Ts) && (t<Tr)) {
				stf_out[i] = r3(t[i], Tr, Ts);
			}

			// tr < t < tr + ts
			else if ((t>=Tr)&&(t<Tr+Ts)) {
				stf_out[i] = r4(t[i], Tr, Ts);
			}

			else if ((t>Tr+Ts)&&(t<Tr+2*Ts)) {
				stf_out[i] = r5(t[i], Tr, Ts);
			}
			// ts + ts < t < tr + 2*ts

		// set of second ranges

			// 0 < t < ts

			// ts < t < tr

			// tr < t < 2*ts

			// 2*ts < t < ts + tr

			// ts + tr < t < tr + 2*ts
		}	
	}
}

/* helping functions to build the tinti function, r stands for range. */
float r1(float t, float Tr) {
	float p1;

	p1 = (Pi*(4*t - Tr)*Tr + 8*pow(t,1.5)*sqrt(-t + Tr) + 4*Tr*sqrt(t*(-t + Tr)) - 
     2*(4*t - Tr)*Tr*atan((-2*t + Tr)/(2.*sqrt(t*(-t + Tr)))))/16.;

	return p1;
}

float r2(float t, float Tr, float Ts) {
	float p1;
	float p2;

	p1 = (Tr*sqrt(t*(-t + Tr)) + 2*sqrt(pow(t,3)*(-t + Tr)) - t*Tr*sqrt((-t + Tr + Ts)/(t - Ts)) + 
      Tr*Ts*sqrt((-t + Tr + Ts)/(t - Ts)) - 2*t*sqrt((t - Ts)*(-t + Tr + Ts)) - 
      2*Ts*sqrt((t - Ts)*(-t + Tr + Ts)) + (4*t - Tr)*Tr*asin(sqrt(t/Tr)) + 
      Tr*(-4*t + Tr)*asin(sqrt((t - Ts)/Tr)))/4.;

	p2 = (-((2*t + Tr - 6*Ts)*sqrt((t - Ts)*(-t + Tr + Ts))) + 
      Tr*(-4*t + Tr + 8*Ts)*asin(sqrt((t - Ts)/Tr)))/4.;

	return p1+p2;
}

float r3(float t, float Tr, float Ts) {
	float p1;
	float p2;

	p1 = (Tr*sqrt(t*(-t + Tr)) + 2*sqrt(pow(t,3)*(-t + Tr)) - t*Tr*sqrt((-t + Tr + Ts)/(t - Ts)) + 
      Tr*Ts*sqrt((-t + Tr + Ts)/(t - Ts)) - 2*t*sqrt((t - Ts)*(-t + Tr + Ts)) - 
      2*Ts*sqrt((t - Ts)*(-t + Tr + Ts)) + (4*t - Tr)*Tr*asin(sqrt(t/Tr)) + 
      Tr*(-4*t + Tr)*asin(sqrt((t - Ts)/Tr)))/4.;

	p2 = (2*(-2*t*sqrt((t - Ts)*(-t + Tr + Ts)) - Tr*sqrt((t - Ts)*(-t + Tr + Ts)) + 
      6*Ts*sqrt((t - Ts)*(-t + Tr + Ts)) + 2*t*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)) + 
      Tr*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)) - 4*Ts*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts))) - 
      Tr*(-4*t + Tr + 8*Ts)*atan((-2*t + Tr + 2*Ts)/(2.*sqrt((t - Ts)*(-t + Tr + Ts)))) + 
      Tr*(-4*t + Tr + 8*Ts)*atan((-2*t + Tr + 4*Ts)/(2.*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)))))/8.;

	return p1+p2;
}

float r4(float t, float Tr, float Ts) {
	float p1;
	
	p1 = t*Tr*acos(sqrt((t - Ts)/Tr)) + (-(sqrt((t - Ts)*(-t + Tr + Ts))*(2*t + Tr + 2*Ts)) - 
      pow(Tr,2)*acos(1.0/sqrt(Tr/(t - Ts))))/4.;

    p2 = (2*(-2*t*sqrt((t - Ts)*(-t + Tr + Ts)) - Tr*sqrt((t - Ts)*(-t + Tr + Ts)) + 
         6*Ts*sqrt((t - Ts)*(-t + Tr + Ts)) + 2*t*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)) + 
         Tr*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)) - 4*Ts*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts))) - 
      Tr*(-4*t + Tr + 8*Ts)*atan((-2*t + Tr + 2*Ts)/(2.*sqrt((t - Ts)*(-t + Tr + Ts)))) + 
      Tr*(-4*t + Tr + 8*Ts)*atan((-2*t + Tr + 4*Ts)/(2.*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)))))/8.;

    return p1+p2;
}

float r5(float t, float Tr, float Ts) {
	float p1;

	p1 = (4*(2*t + Tr - 4*Ts)*sqrt(-((t - 2*Ts)*(t - Tr - 2*Ts))) + Pi*Tr*(-4*t + Tr + 8*Ts) + 
     	  2*Tr*(-4*t + Tr + 8*Ts)*atan((-2*t + Tr + 4*Ts)/(2.*sqrt((t - 2*Ts)*(-t + Tr + 2*Ts)))))/16.;

   return p1;
}

/* used for different integration range when Ts < Tr < 2 Ts */
float sr3(float* ts, float Tr, float Ts) {
	float p1, p2;

	p1 = t*Tr*acos(sqrt((t - Ts)/Tr)) + (-(sqrt((t - Ts)*(-t + Tr + Ts))*(2*t + Tr + 2*Ts)) - 
      pow(Tr,2)*acos(1.0/sqrt(Tr/(t - Ts))))/4.;

	p2 = (-4*(2*t + Tr - 6*Ts)*sqrt((t - Ts)*(-t + Tr + Ts)) + Pi*Tr*(-4*t + Tr + 8*Ts) - 
      2*Tr*(-4*t + Tr + 8*Ts)*atan((-2*t + Tr + 2*Ts)/(2.*sqrt((t - Ts)*(-t + Tr + Ts)))))/16.;

	return p1+p2;
}