/* This program creates a kinematic source from results of an AWP dynamic rupture
   simulation, both for a vertical non-planar fault.

    William Savran, wsavran@ucsd.edu

   MPI-IO is used both for reading the source time function and writing
   the moment rate file.
*/

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include "sord_mpio.h"
#include "utils.h"
#include "xapiir.h"
#include "spline.h"

/* TODO (later): Write de-staggering routine as option
 *   (1) implement parameter file */


/* computed scalar moment given moment rate tensor */
float compute_moment(int nst, float dt, float *xx, float *yy, float *zz, float *xz, float *yz, float *xy);

int main (int argc, char*argv[]) {
    /* modify these parameters */
    int nchunks = 1; /* also equal to the number of i/o calls */
    int nx = 2601, ny = 801; 
    float dt_sord = 0.002, dt_awp = 0.00175;
    float rupture_time = 20.0010;
    int nt_sord = 10001;
    int dx_sord = 25;
    float mean_faultn_coord = 20000.0;
    int x_start = 1600, y_start = 1600, z_start = 1;/* node locations */
    char *nhat1_file = "./out/nhat1";
    char *nhat2_file = "./out/nhat2";
    char *nhat3_file = "./out/nhat3";
    char *momrate_file = "./striped/momrate";
    char *mu_file = "./out/mu_two_step";
    char *ycomp_file = "./out/x3o";
    char *sv1_file = "./out/sv1";
    char *sv2_file = "./out/sv2";
    char *sv3_file = "./out/sv3";

    int debug = 0;
    int filter = 1;
    int interp = 1; /* 0: false, 1:true */
    int partition = 0; /* 0: do not write out tpsrc_tmp, 1: write tpsrc_tmp */

    /*parameters for xapiir*/
    int iord=4, npas=2;
    float trbndw=0., a=0.; /* chebyshev parameters */
    char *aproto="BU";
    char *ftype="LP";
    float hp=0., lp=10.00; /*lp=1.0*/
    
    /*variables for splint*/
    float *y2_sv1, *y2_sv2, *y2_sv3;
    float yp1 = 0.0; float ypn = 0.0;
    
    /* Don't modify anything below here! */
    float *buf_sv1, *buf_sv2, *buf_sv3;
    float **sv1, **sv2, **sv3;
    float *buf_nhat1, *buf_nhat2, *buf_nhat3;
    float *buf_mu, *buf_ycomp;
    int nt_awp;
    float *t_awp, *t_sord;
    int csize;
    float sv1_awp, sv2_awp, sv3_awp;
    int rank, nprocs;
    int s0, xi, yi;
    int k, n, l, i;
    float **xx, **yy, **zz, **xz, **yz, **xy;
    int *xil, *yil, *zil;
    float area;
    int temp, nst;
    MPI_Offset off;
    float global_moment;
    float moment, temp_moment;
    
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
    buf_mu = (float*)calloc(csize, sizeof(float));
    buf_ycomp = (float*)calloc(csize, sizeof(float));
    xil = (int*)calloc(csize, sizeof(int));
    yil = (int*)calloc(csize, sizeof(int));
    zil = (int*)calloc(csize, sizeof(int));
    
    /* 2D arrays for time-series */
    sv1 = (float**)calloc(csize, sizeof(float*));
    sv2 = (float**)calloc(csize, sizeof(float*));
    sv3 = (float**)calloc(csize, sizeof(float*));
    
    for (l=0; l<csize; l++) {
        sv1[l] = (float*) calloc(nt_sord, sizeof(float));
        sv2[l] = (float*) calloc(nt_sord, sizeof(float));
        sv3[l] = (float*) calloc(nt_sord, sizeof(float));
    }
    
    /* if we are interpolating, compute time vectors,
     * check if t_sord is equal to nt_sord */
    if (interp) {
        t_sord = arange(0.0, rupture_time, dt_sord, &temp);
        t_awp = arange(0.0, rupture_time, dt_awp, &nt_awp);
        nst = nt_awp;
        if (nt_sord != temp) {
            fprintf(stderr,"%d", temp);
            fprintf(stderr,"interpolation time-vector length not equal\n");
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Finalize();
            exit(1); 
        }
        fprintf(stderr, "nst: %i\n", nst); 
        /* 2nd derivatives of slip rates for cubic splines */
        y2_sv1 = (float*)calloc(nt_sord, sizeof(float));
        y2_sv2 = (float*)calloc(nt_sord, sizeof(float));
        y2_sv3 = (float*)calloc(nt_sord, sizeof(float));
    }
    else {
        nst = nt_sord;
        if (rank==0) fprintf(stderr,"nst: %d\n", nst);
    }
    
    // allocate arrays
    xx = (float**)calloc(csize, sizeof(float*));
    yy = (float**)calloc(csize, sizeof(float*));
    zz = (float**)calloc(csize, sizeof(float*));
    xz = (float**)calloc(csize, sizeof(float*));
    yz = (float**)calloc(csize, sizeof(float*));
    xy = (float**)calloc(csize, sizeof(float*));
    
    for (i=0; i<csize; i++) {
        xx[i] = (float*)calloc(nst, sizeof(float));
        yy[i] = (float*)calloc(nst, sizeof(float));
        zz[i] = (float*)calloc(nst, sizeof(float));
        xz[i] = (float*)calloc(nst, sizeof(float));
        yz[i] = (float*)calloc(nst, sizeof(float));
        xy[i] = (float*)calloc(nst, sizeof(float));
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
        read_fault_params(mu_file, off, csize, buf_mu);
        read_fault_params(ycomp_file, off, csize, buf_ycomp);

        /* organize into 2d arrays for easy filtering */
        for (n=0; n < nt_sord; n++) {
            for (l=0; l < csize; l++) {
                sv1[l][n]=buf_sv1[n*csize+l];
                sv2[l][n]=buf_sv2[n*csize+l];
                sv3[l][n]=buf_sv3[n*csize+l];
            }
        }
        for (l=0; l<csize; l++) {
            /* filter slip-rate functions, 
             * maybe filter after? */
            if (filter) {
                xapiir_(sv1[l], &nt_sord, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt_sord, &npas);
                xapiir_(sv2[l], &nt_sord, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt_sord, &npas);
                xapiir_(sv3[l], &nt_sord, aproto, &trbndw, &a, &iord, ftype, &hp, &lp, &dt_sord, &npas);
            }
            
            /* either interpolate or don't... */
            area = dx_sord * dx_sord;
            // mu = 3464.0 * 3464.0 * 3000.0; 
            if (interp) {
                /* just some output */
                if (rank==0 && l == 0) fprintf(stderr,"interpolating time-series. . .\n");
                if (rank==0 && l == 0) fprintf(stderr,"old dt: %f\nnew dt: %f\nnt: %d\n", dt_sord, dt_awp, nst);

                /* interpolate slip rates to same dt as awp simulation */
                spline(t_sord, sv1[l], nt_sord, yp1, ypn, y2_sv1);
                spline(t_sord, sv2[l], nt_sord, yp1, ypn, y2_sv2);
                spline(t_sord, sv3[l], nt_sord, yp1, ypn, y2_sv3);

                for (n=0; n<nt_awp; n++) {
                    splint(t_sord, sv1[l], y2_sv1, nt_sord, t_awp[n], &sv1_awp);
                    splint(t_sord, sv2[l], y2_sv2, nt_sord, t_awp[n], &sv2_awp);
                    splint(t_sord, sv3[l], y2_sv3, nt_sord, t_awp[n], &sv3_awp);
                    
                    /* mom = outer product of slip, nhat (potency) * area * shear modulus.
                       see Shi and Day, 2013 for information on SORD coordinate system 
                    xx[l][n] = buf_mu[l] * area * (sv1[l][n] * buf_nhat1[l] + sv1[l][n] * buf_nhat1[l]);  
                    yy[l][n] = buf_mu[l] * area * (sv3[l][n] * buf_nhat3[l] + sv3[l][n] * buf_nhat3[l]); 
                    zz[l][n] = buf_mu[l] * area * (sv2[l][n] * buf_nhat2[l] + sv2[l][n] * buf_nhat2[l]); 
                    xz[l][n] = buf_mu[l] * area * (sv1[l][n] * -buf_nhat2[l] - sv2[l][n] * buf_nhat1[l]); 
                    yz[l][n] = buf_mu[l] * area * (sv3[l][n] * -buf_nhat2[l] - sv2[l][n] * buf_nhat3[l]); 
                    xy[l][n] = buf_mu[l] * area * (sv1[l][n] * buf_nhat3[l] + sv3[l][n] * buf_nhat1[l]); */
                    xx[l][n] = buf_mu[l] * area * (sv1_awp * buf_nhat1[l] + sv1_awp * buf_nhat1[l]);  
                    yy[l][n] = buf_mu[l] * area * (sv3_awp * buf_nhat3[l] + sv3_awp * buf_nhat3[l]); 
                    zz[l][n] = buf_mu[l] * area * (sv2_awp * buf_nhat2[l] + sv2_awp * buf_nhat2[l]); 
                    xz[l][n] = buf_mu[l] * area * (sv1_awp * -buf_nhat2[l] - sv2_awp * buf_nhat1[l]); 
                    yz[l][n] = buf_mu[l] * area * (sv3_awp * -buf_nhat2[l] - sv2_awp * buf_nhat3[l]); 
                    xy[l][n] = buf_mu[l] * area * (sv1_awp * buf_nhat3[l] + sv3_awp * buf_nhat1[l]);
                }                                

                
                /* compute moment */
                temp_moment = compute_moment(nst, dt_awp, xx[l], yy[l], zz[l], xz[l], yz[l], xy[l]);
                moment += temp_moment;
            }
            /* no interpolation */
            else {
                for (n=0; n<nt_sord; n++) {
                    // mom = outer product of potency and material properties 
                    // These coordinates are in SORD format, ie x1 = x; x2 = z; x3 = y
                    xx[l][n] = buf_mu[l] * area * (sv1[l][n] * buf_nhat1[l] + sv1[l][n] * buf_nhat1[l]);  
                    yy[l][n] = buf_mu[l] * area * (sv3[l][n] * buf_nhat3[l] + sv3[l][n] * buf_nhat3[l]); 
                    zz[l][n] = buf_mu[l] * area * (sv2[l][n] * buf_nhat2[l] + sv2[l][n] * buf_nhat2[l]); 
                    xz[l][n] = buf_mu[l] * area * (sv1[l][n] * -buf_nhat2[l] - sv2[l][n] * buf_nhat1[l]); 
                    yz[l][n] = buf_mu[l] * area * (sv3[l][n] * -buf_nhat2[l] - sv2[l][n] * buf_nhat3[l]); 
                    xy[l][n] = buf_mu[l] * area * (sv1[l][n] * buf_nhat3[l] + sv3[l][n] * buf_nhat1[l]); 
                }     
                /* compute moment */
                temp_moment = compute_moment(nt_sord, dt_sord, xx[l], yy[l], zz[l], xz[l], yz[l], xy[l]);
                moment += temp_moment;
            }
            /* determine subfault location assuming fault is striking along x component */
            xil[l] = x_start + (s0+l) % nx; 
            yil[l] = y_start + rint((buf_ycomp[l] - mean_faultn_coord) / dx_sord);
            zil[l] = z_start + (s0+l) / nx;
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
        free(xx[i]);
        free(yy[i]);
        free(zz[i]);
        free(xz[i]);
        free(yz[i]);
        free(xy[i]);
        free(sv1[i]);
        free(sv2[i]);
        free(sv3[i]);
    }
    free(xx);
    free(yy);
    free(zz);
    free(xz);
    free(yz);
    free(xy);
    free(sv1);
    free(sv2);
    free(sv3);
    free(buf_nhat1);
    free(buf_nhat2);
    free(buf_nhat3);
    free(buf_sv1);
    free(buf_sv2);
    free(buf_sv3);
    free(buf_ycomp);
    free(buf_mu);
    free(xil);
    free(yil);
    free(zil);
    if (interp) {
        free(y2_sv1);
        free(y2_sv2);
        free(y2_sv3);
    } 

    /* finalize mpi */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

     
float compute_moment(int nst, float dt, float *xx, float *yy, float *zz, float *xz, float *yz, float *xy) {
    float mxx, myy, mzz, mxz, myz, mxy;
    float mij_squared;
    float m0;
    float prefactor;

    /* integrate moment-rate to get moment */
    mxx = sum(xx, nst)*dt;
    myy = sum(yy, nst)*dt;
    mzz = sum(zz, nst)*dt;
    mxz = sum(xz, nst)*dt;
    myz = sum(yz, nst)*dt;
    mxy = sum(xy, nst)*dt;

    /* square moment rate tensor */
    prefactor = 1.0f / sqrt(2);
    mij_squared = mxx*mxx + myy*myy + mzz*mzz + 2*mxz*mxz + 2*myz*myz + 2*mxy*mxy;
    m0 = prefactor * sqrt(mij_squared);
    return m0;
}
