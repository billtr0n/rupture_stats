/* This program performs fitting of slip-rate functions for statistical analysis of source fields.

    William Savran, wsavran@ucsd.edu
    Date: 04.26.2016

   MPI-IO is used both for reading the source time function and writing
   the moment rate file.
*/

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include "sord_mpio.h"
#include "utils.h"
#include "brute.h"
#include "xapiir.h"
#include "spline.h"


int main (int argc, char*argv[]) {
    /* modify these parameters */
    int nchunks = 1; /* also equal to the number of i/o calls */
    int nx = 2601, ny = 801; int nt=10001;
    float dt = 0.002;
    int dx = 25;
    char *nhat1_file = "./out/nhat1";
    char *nhat2_file = "./out/nhat2";
    char *nhat3_file = "./out/nhat3";
    char *sv1_file = "./out/sv1";
    char *sv2_file = "./out/sv2";
    char *sv3_file = "./out/sv3";
    
    /* Don't modify anything below here. */
    float *buf_sv1, *buf_sv2, *buf_sv3;
    float **sv1, **sv2, **sv3; **svm;
    float *buf_nhat1, *buf_nhat2, *buf_nhat3;
    float *t_sord;

    int csize;
    int rank, nprocs;
    int s0, xi, yi;
    int k, n, l, i;
    
    float *buf_tp; *buf_tr;
   
    MPI_Offset off;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    

    /* number of points to process per read operation */
    if (((nx*ny) % (nprocs * nchunks)) != 0) {
        fprintf(stdout, "number of points not divisible by number of cpus * buffer size\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(1);
    } 
    
    /* number of time-series read per mpi-io call */
    csize = nx*ny / nprocs / nchunks;
    
    /* allocate arrays */
    buf_sv1 = (float*)calloc(csize*nt_sord, sizeof(float));
    buf_sv2 = (float*)calloc(csize*nt_sord, sizeof(float));
    buf_sv3 = (float*)calloc(csize*nt_sord, sizeof(float));
    buf_nhat1 = (float*)calloc(csize, sizeof(float));    
    buf_nhat2 = (float*)calloc(csize, sizeof(float));
    buf_nhat3 = (float*)calloc(csize, sizeof(float));
    buf_tp = (float*)calloc(csize*sizeof(float));
    buf_tr = (float*)calloc(csize*sizeof(float));
    
    /* 2D arrays for time-series */
    sv1 = (float**)calloc(csize, sizeof(float*));
    sv2 = (float**)calloc(csize, sizeof(float*));
    sv3 = (float**)calloc(csize, sizeof(float*));
    svm = (float**)calloc(csize, sizeof(float*));
    
    for (l=0; l<csize; l++) {
        sv1[l] = (float*) calloc(nt_sord, sizeof(float));
        sv2[l] = (float*) calloc(nt_sord, sizeof(float));
        sv3[l] = (float*) calloc(nt_sord, sizeof(float));
        svm[l] = (float*) calloc(nt_sord, sizeof(float));
    }

    /* main working loop */
    MPI_Barrier(MPI_COMM_WORLD);
    moment = 0;
    for (k = 0; k < nchunks; k++) {
        s0 = rank*nchunks*csize + k*csize;
        off = (MPI_Offset) s0 * sizeof(float);
        yi = s0 / nx;
        xi = s0 % nx;   
        if (debug == 1) {
            // output some parameters
            fprintf(stderr, "rank: %i\n", rank);
            fprintf(stderr, "s0: %i\n", s0);
            fprintf(stderr, "off: %i\n", off);
            fprintf(stderr, "yi: %i\n", yi);
            fprintf(stderr, "xi: %i\n\n", xi);
            
        }
                
        /* read dynamic rupture parameters */ 
        if (rank == 0) fprintf(stderr, "Reading time-series...\n");
        read_time_series(sv1_file, nx, ny, xi, yi, nt_sord, csize, buf_sv1);
        read_time_series(sv2_file, nx, ny, xi, yi, nt_sord, csize, buf_sv2);
        read_time_series(sv3_file, nx, ny, xi, yi, nt_sord, csize, buf_sv3);

        /* these routines might not need to be mpi/io files are quite small
           might be better to read in at rank zero process and broadcast */
        read_fault_params(nhat1_file, off, csize, buf_nhat1);
        read_fault_params(nhat2_file, off, csize, buf_nhat2);
        read_fault_params(nhat3_file, off, csize, buf_nhat3);

        /* organize into 2d arrays for easy filtering */
        for (n=0; n < nt_sord; n++) {
            for (l=0; l < csize; l++) {
                sv1[l][n]=buf_sv1[n*csize+l];
                sv2[l][n]=buf_sv2[n*csize+l];
                sv3[l][n]=buf_sv3[n*csize+l];
                svm[l][n]=sqrt(powf(sv1[l][n],2)+powf(sv2[l][n],2)+powf(sv3[l][n],2));
            }
        }

        /* perform evalulations of slip-rate files */
        for (l=0; l<csize; l++) {

        	    
        } /*end csize loop */
        /* write out moment rates. this implementation writes out once per chunk. keep in
           mind for memory. */
        if (rank==0) fprintf(stderr,"accumulated moment: %f\n", moment);
        if (rank==0) fprintf(stderr,"writing moment rates. . . %i chunks out of %i \n", k+1, nchunks);
        write_momrate(momrate_file, nst, nchunks, rank, csize, k, xil, yil, zil, xx, yy, zz, xz, yz, xy);
        // write locations to tpsrc_tmp for source partitioner
        if (partition==1) {
            if (rank==0) fprintf(stderr, "writing subfault locations...\n");
            write_locations("tpsrc_tmp", off, csize, xil, yil, zil);
        }
    } /* end main loop */

    /* compute moment from all processes */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&moment, &global_moment, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
    if (rank == 0) {
        fprintf(stderr, "total moment: %e\n", global_moment);
        fprintf(stderr, "mw: %e\n", (2.0/3.0)*(log10(global_moment)-9.1));
    }
    /* free buffers */
    for (i=0;i<csize;i++) {
        free(sv1[i]);
        free(sv2[i]);
        free(sv3[i]);
        free(svm[i]);
    }
    free(sv1);
    free(sv2);
    free(sv3);
    free(svm);
    free(buf_nhat1);
    free(buf_nhat2);
    free(buf_nhat3);
    free(buf_sv1);
    free(buf_sv2);
    free(buf_sv3);

    /* finalize mpi */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}