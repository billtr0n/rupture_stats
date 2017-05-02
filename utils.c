#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

float* linspace(float start, float end, int length) {
    float* buf;
    float dx;
    int i;
    
    dx = (end - start) / (length-1);
    buf = (float*)malloc(length*sizeof(float));
    for (i=0; i<length; i++) { 
        buf[i] = start + dx*i;
    }
    return buf;
}

float* arange(float start, float end, float dx, int *n) {
    float* buf;
    int length;
    int i;
    
    length = ceil( (end - start) / dx);
    buf = (float*)malloc(length*sizeof(float));
    for (i=0; i<length; i++) {
        buf[i] = start + dx*i;
    }
    *n = length;
    return buf;

}

float* zeros_f(int n) {
    float* buf;
    int i;

    buf = (float*)malloc(n*sizeof(float));
    for (i=0; i<n; i++) {
        buf[i] = 0.0;
    }

    return buf;
}

int* zeros_i(int n) {
    int* buf;
    int i;

    buf = (int*)malloc(n*sizeof(int));
    for (i=0; i<n; i++) {
        buf[i] = 0;
    }

    return buf;
}

float sum(float *buf, int n) {
    float total;
    int i;
    total = 0;
    for (i=0;i<n;i++) {
        total += buf[i];
    }

    return total;
}

float** cartesian_product( float **A, int *np, int *ncom, int ndim ) {
    int i, j;
    int n = 1;
    bool DEBUG = false;
    int *c;
    float **B;

    // compute ncom
    for (i=0; i<ndim; i++) {
        n *= np[i];
    }
    *ncom = n;

    // allocate count array
    c = zeros_i( ndim );

    // allocate output array
    B = (float**) malloc( sizeof(float*) * n );
    for (i=0; i<n; i++) {
        B[i] = (float*) malloc(sizeof(float)*ndim);
    }

    // prepare combinations
    if (DEBUG) {
        printf("\n");
        printf("Combinations\n");
        printf("============\n");
    }

    for (i=0; i<n; i++) {
        for (j=0; j<ndim; j++) {
            if (c[j]==np[j]) {
                c[j+1]++;
                c[j]=0;
            }
            B[i][j] = A[j][c[j]];
        }
        if (DEBUG) {
            printf("{%d, %d, %d} -> {%f, %f %f}\n", c[0], c[1], c[2], B[i][0], B[i][1], B[i][2]);
        }
        
        c[0]++;
    }
        
    // all done!
    return B;

}


bool compare_f( float a, float b, float eps ) {
    return fabs(a-b) < eps;       
}

void dealloc2d_f( float** buf, int n ) {
    int i;

    for (i=0; i<n; i++) {
        free(buf[i]);
    }
    free(buf);
}

void dealloc2d_i( int** buf, int n ) {
    int i;

    for (i=0; i<n; i++) {
        free(buf[i]);
    }
    free(buf);
}

float **outer( float *v1, float *v2, int n ) {
    int i, j;
    float **m;

    m = (float**) malloc(sizeof(float*) * n);
    for (i=0; i<n; i++) {
        m[i] = (float*) malloc(sizeof(float) * n);
    }

    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            m[i][j] = v1[i] * v2[j];
        }
    }

    return m;

}

float rastrigin( float x, float y ) {
    float f;
    float p1,p2,p3;
    int i, j;
    int d = 2;
    float pi = 3.14159265359;

    p1 = 10*d;
    p2 = x*x + y*y;
    p3 = 10 * (cos(2*pi*x) + cos(2*pi*y));
    
    f = p1 + p2 - p3;
    return f;
}

float ackley( float x, float y ) {
    float p1, p2, f;
    float pi = 3.14159265359;

    p1 = 20*exp(-0.2 * sqrt(0.5*(x*x+y*y)));
    p2 = exp(0.5*(cos(2*pi*x) + cos(2*pi*y)));
    f = -p1 - p2 + exp(1) + 20;
    return f;
}


/*
int main() {
    float *buf;
    float sum1;
    int i;
    buf = (float*)malloc(10*sizeof(float));

    for (i=0; i<10; i++) {
        buf[i] = i*1e8;
    }
    sum1 = sum(buf, 10);

    fprintf(stdout, "%f\n", sum1);

    return 0;
}
    
int main() {
    float start = 0;
    float end = 10;
    float dx = 1;
    int length = 11;
    float *buf1, *buf2;
    int i;
    int temp;
    
    buf1 = linspace(start, end, length);
    buf2 = arange(start, end+0.1, dx, &temp);

    
    for (i=0; i<length; i++) {
        printf("%f %f\n", buf1[i], buf2[i]);
    } 
    
    free(buf1);
    free(buf2);
    return 0;
}
*/
