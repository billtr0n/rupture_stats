#include <stdio.h>
#include <math.h>
#include "utils.h"

#define FTOL 1e-7f
#define PI 3.14159265359
#define RAD2DEG 57.2957795131

float get_dip(float nhat1, float nhat2, float nhat3) {

 	float dip;
 	int i, j;
 	float nproj[3], n[3];
 	float scaling;
 	float arg;
 	float theta;
 	
 	// assumes that free surface is at x1-x3 plane
    nproj[0] = nhat1;
    nproj[1] = 0;
    nproj[2] = nhat3;
    to_unit_vector(nproj, 3);


    // assembles normal vector from fault
    n[0] = nhat1;
    n[1] = nhat2;
    n[2] = nhat3;
    to_unit_vector(n, 3);

    // argument of dot product
    arg = ( n[0]*nproj[0] + n[2]*nproj[2] );

    theta = acos(arg) * RAD2DEG;
   
    if (nhat2 > 0) {
        dip = 90 + theta;
    }
    else if (nhat2 <= 0) {
    	dip = 90 - theta;
    }
    
 	return dip;
}

float get_strike(float nhat1, float nhat3) {
	float nproj[3];
	float surf_normal[3];
	float strike;
    float theta;

    // implicit surface in x1-x3 plane
    nproj[0] = nhat1;
    nproj[1] = 0.0f;
    nproj[2] = nhat3;
    to_unit_vector(nproj, 3);
    // printf("nproj[0]=%f nproj[1]=%f nproj[2]=%f\n", nproj[0], nproj[1], nproj[2]);

    // implicit multiplication of nproj[0] with x1_hat
    // x1_hat = [1, 0, 0]
    theta = acos( nproj[0] ) * RAD2DEG;
    // printf("theta=%f\n", theta);
    
    // convention so that if the fault is striking toward positive x3 the strike is positive
    strike = theta - 90.0;
    return strike;
}

float get_rake(float nhat1, float nhat2, float nhat3, float su1, float su2, float su3) {
    int i;
    float s;
    float slip[3], ns[3];
    float arg;
    float scaling;
    float rake;

    // stored in seperate files. avoiding pointers if necessary.
    slip[0] = su1;
    slip[1] = su2;
    slip[2] = su3;

    // check if fault ruptured
    s = norm( slip, 3 );
    if ((s - FTOL) <= 0) {
        return -99999.0;
    }

    // nhat x [0 1 0]
    // (nhat[2]) - (nhat[0])
    ns[0] = nhat3;
    ns[1] = 0;
    ns[2] = -nhat1;
    to_unit_vector(ns, 3);
    to_unit_vector(slip, 3);

    // cross product
    arg = slip[0]*ns[0]+slip[1]*ns[1]+slip[2]*ns[2];

    //printf("arg=%f scaling=%f\n", arg, scaling);
    rake = acos( arg ) * RAD2DEG;
    return rake;
}

/*
int main() {
    float nhat1, nhat2, nhat3;
    float su1, su2, su3;
    float strike, dip, rake;

    nhat1=-0.057126;
    nhat2=0.082879;
    nhat3=0.997446;

    nhat1 = 0;
    nhat2 = 0;
    nhat3 = 1;

    //nhat2=-0.221027;
    //nhat1=-0.022632;
    //nhat3=0.975005;

    //su1=0.998367;
    //su2=-0.001967;
    //su3=0.057093;
    
    su1=1;
    su2=0;
    su3=0;

    strike = get_strike(nhat1, nhat3);
    dip = get_dip(nhat1, nhat2, nhat3);
    rake = get_rake(nhat1, nhat2, nhat3, su1, su2, su3);
    
    //printf("Test 1: Strike and Slip in same direction.\n\tExpect Strike = 0.0 Dip = 90.0 Rake = 0.0\n");
    printf("nhat = [%f %f %f]\nslip = [%f %f %f]\n", nhat1, nhat2, nhat3, su1, su2, su3);
    printf("strike=%f\ndip=%f\nrake=%f\n\n", strike, dip, rake);

    strike = get_strike(nhat1,nhat3);
    dip = get_dip(nhat1, nhat2, nhat3);
    rake = get_rake(nhat1, nhat2, nhat3, su1, su2, su3);
    //printf("Test 2: Strike along x1 and Slip = [1, 0, 1].\n\tExpect Strike = 0.0 Dip = 90.0 Rake = 45.0\n");
    printf("nhat = [%f %f %f]\nslip = [%f %f %f]\n", nhat1, nhat2, nhat3, su1, su2, su3);
    printf("strike=%f\ndip=%f\nrake=%f\n\n", strike, dip, rake);
    
    return 0;
}
*/
