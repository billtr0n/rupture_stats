#include <stdio.h>
#include "../utils.h"
#include <stdlib.h>


int main() {
	float **buf = NULL;
	int i;
	int n = 100;

    if (true == NULL) {
        printf("pointer is NULL\n");
    }
	buf = (float**) malloc( sizeof(float*) * n);
	for (i=0; i<n; i++) {
		buf[i] = (float*) malloc( sizeof(float) * 10 );
	}
	dealloc2d_f( buf, n );


	return 0;
}
