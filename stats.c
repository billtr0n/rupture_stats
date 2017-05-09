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
#include "stf.h"

#define ndim 1


int main (int argc, char*argv[]) {
    //////////////////////////////////////////////////////////////////////////
    // USER PARAMETERS
    //////////////////////////////////////////////////////////////////////////
    int nchunks = 1; /* also equal to the number of i/o calls */
    int nx = 2601, ny = 801; int nt=10001;
    float dt = 0.002;
    int dx = 25;
    float pl[ndim], ph[ndim], pd[ndim];
    
    // define parameters for estimation
    pl[0] = dt;
    ph[0] = Tr/2;
    pd[0] = dt;

    // input file names
    char *nhat1_file = "./out/nhat1";
    char *nhat2_file = "./out/nhat2";
    char *nhat3_file = "./out/nhat3";
    char *sv1_file = "./out/sv1";
    char *sv2_file = "./out/sv2";
    char *sv3_file = "./out/sv3";

    // output file names
    char *tp_file = "./out/tp";
    char *tr_file = "./out/tr";
    char *slip_file = "./out/slip";
    char *slip_kin_file = "./out/slip_kin";
    char *psv_file = "./out/psv";
    char *psv_kin_file = "./out/psv_kin";
    char *t0_file = "./out/t0";
    char *tarr_file = "./out/tarr";
    char *tarr_kin_file = "./out/tarr_kin";
    char *dip_file = "./out/dip";
    char *rake_file = "./out/rake";
    char *strike_file = "./out/strike";
    
    //////////////////////////////////////////////////////////////////////////
    // Don't modify anything below here. 
    //////////////////////////////////////////////////////////////////////////
    obj_fun_ptr func = &objective_function;

    float *buf_sv1, *buf_sv2, *buf_sv3;
    float **sv1, **sv2, **sv3; **svm;
    float *buf_nhat1, *buf_nhat2, *buf_nhat3;

    float *buf_tp; *buf_tr;
    float *buf_slip_kin, *buf_psv, *buf_psv_kin, *buf_t0, *buf_tarr, *buf_slip;
    float *buf_strike, *buf_dip, *buf_rake;
    float *t_sord;

    float su1, su2, su3;
    float *su;
    float kin_slip;
    float tarr;
    float tarr_kin;
    float **comb, **par;
    float *t, *stf;

    int tarr_ind, tarr_kin_ind;
    int tmp_n;
    int n_comb;
    int *n_par;

    int csize;
    int rank, nprocs;
    int s0, xi, yi;
    int k, n, l, i;
    
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
    buf_sv1 = calloc(csize*nt_sord, sizeof(float));
    buf_sv2 = calloc(csize*nt_sord, sizeof(float));
    buf_sv3 = calloc(csize*nt_sord, sizeof(float));
    buf_nhat1 = calloc(csize, sizeof(float));    
    buf_nhat2 = calloc(csize, sizeof(float));
    buf_nhat3 = calloc(csize, sizeof(float));

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

    // calculated arrays
    buf_tp = calloc(csize*sizeof(float));
    buf_tr = calloc(csize*sizeof(float));
    buf_strike = calloc(csize*sizeof(float));
    buf_dip = calloc(csize*sizeof(float));
    buf_rake = calloc(csize*sizeof(float));
    buf_slip = calloc(csize*sizeof(float));
    buf_slip_kin = calloc(csize*sizeof(float));
    buf_psv = calloc(csize*sizeof(float));
    buf_psv_kin = calloc(csize*sizeof(float));
    buf_t0 = calloc(csize*sizeof(float));
    buf_tarr = calloc(csize*sizeof(float));
    buf_tarr_kin = calloc(csize*sizeof(float));

    n_par = malloc(sizeof(int) * ndim);
    par = malloc(sizeof(float*) * ndim);

    // create possibly jagged array of parameters       
    for (i=0; i<ndim; i++)  {
        par[i] = arange(pl[i], ph[i], pd[i], &n_par[i]);
    }

    // generate all possible combinations of parameters
    comb = cartesian_product( par, n_par, &n_comb, ndim );

    // to store full_output of cost function
    full = malloc( sizeof(float) * n_comb);

    // time vector
    t = arange(0, dt*nt, dt, &tmp_n);
    if (nt != tmp_n) {
        printf("ERROR: sim nt and stf nt not equivalent. aborting!\n");
        MPI_Finalize();
        exit(-1);
    }

    /* main working loop */
    MPI_Barrier(MPI_COMM_WORLD);
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
                sv1[l][n] = buf_sv1[n*csize+l];
                sv2[l][n] = buf_sv2[n*csize+l];
                sv3[l][n] = buf_sv3[n*csize+l];
                svm[l][n] = sqrt(powf(sv1[l][n],2)+powf(sv2[l][n],2)+powf(sv3[l][n],2));
            }
        }

        // perform evalulations of slip-rate files
        for (l=0; l<csize; l++) {
            
            // compute rupture onset
            t0_ind = find_first( svm[l], nt, ">", 0.001 );
            buf_t0[l] = t0_ind * dt;

            // compute kinematic rupture arrest
            // note: defined as first occurence after t0 where svm drops below 0.001
            tarr_kin_ind = t0_ind + find_first( &svm[l][t0_ind+1], nt-t0_ind-1, "<", 0.001 );
            buf_tarr_kin[l] = tarr_kin_ind * dt;

            // compute dynamic rupture arrest
            // note: defined where 
            tarr_ind = find_last( svm[l], nt, ">", 0.001 );
            buf_tarr[l] = tarr_ind * dt;

            // compute rise time
            buf_tr[l] = tarr_kin - buf_t0[l];
            // printf("t0_ind=%d tarr_kin_ind=%d t0=%f tarr_kin=%f Tr=%f\n", t0_ind, tarr_kin_ind, buf_t0[l], tarr_kin, buf_tr[l]);

            // assemble dynamic slip vector and normalize for angle estimations
            su1 = trapz( sv1[l], nt, dt );
            su2 = trapz( sv2[l], nt, dt );
            su3 = trapz( sv3[l], nt, dt );
            su = unit_vector_3d( su1, su2, su3 );

            // compute kinematic slip
            tmp_n = nt - (t0_ind + ((nt-1) - tarr_ind));
            buf_slip_kin[l] = trapz( &svm[l][t0_ind], tmp_n, dt );

            // compute dynamic slip
            buf_slip[l] = trapz( svm[l], nt, dt );

            // compute peak slip velocity
            buf_psv[l] = maximum( svm[l], nt );
            // printf("k_slip=%f d_slip=%f psv=%f\n", buf_slip_kin[l], buf_slip[l], buf_psv[l]);

            // compute strike, dip, and rake
            buf_strike[l] = get_strike( buf_nhat1[l], buf_nhat3[l] );
            buf_dip[l] = get_dip( buf_nhat1[l], buf_nhat2[l], buf_nhat3[l] );
            buf_rake[l] = get_rake( buf_nhat1[l], buf_nhat2[l], buf_nhat3[l], su[0], su[1], su[2] );

            // set-up extras
            extras[0] = (void*) svm[l];
            extras[1] = (void*) &nt;
            extras[2] = (void*) &dt;
            extras[3] = (void*) &buf_slip_kin[l];
            extras[4] = (void*) &t0_ind;
            extras[5] = (void*) &tarr_kin_ind;
            extras[6] = (void*) &buf_tr[l];

            // printf("n_comb=%d ndim=%d\n", n_comb, ndim);

            // find best fitting peak-time
            min = fmin_brute( func, comb, ndim, n_comb, extras, full );

            // compute tinti function with best fitting parameters
            stf = tinti(t, nt, buf_tr[l], min[0], buf_t0[l], buf_slip_kin[l]);

            // compute kinematic psv
            buf_psv_kin[l] = maximum( stf, nt );            

            // end nicely.
            free(su);

        } /*end csize loop */
            
        /* write out necessary info */
        write_fault_params(tp_file, off, csize, buf_tp, MPI_FLOAT); 
        write_fault_params(tr_file, off, csize, buf_tr, MPI_FLOAT); 
        write_fault_params(slip_file, off, csize, buf_slip, MPI_FLOAT);
        write_fault_params(slip_kin_file, off, csize, buf_slip_kin, MPI_FLOAT); 
        write_fault_params(psv_file, off, csize, buf_psv, MPI_FLOAT); 
        write_fault_params(psv_kin_file, off, csize, buf_psv_kin, MPI_FLOAT); 
        write_fault_params(t0_file, off, csize, buf_t0, MPI_FLOAT); 
        write_fault_params(tarr_file, off, csize, buf_tarr, MPI_FLOAT); 
        write_fault_params(tarr_kin_file, off, csize, buf_tarr_kin, MPI_FLOAT);
        write_fault_params(dip_file, off, csize, buf_dip, MPI_FLOAT);
        write_fault_params(rake_file, off, csize, buf_rake, MPI_FLOAT);
        write_fault_params(strike_file, off, csize, buf_strike, MPI_FLOAT);

    } /* end main loop */

    /* compute moment from all processes */
    MPI_Barrier(MPI_COMM_WORLD);

    /* free buffers */
    dealloc2d_f(sv1, csize);
    dealloc2d_f(sv2, csize);
    dealloc2d_f(sv3, csize);
    dealloc2d_f(svm, csize);
    dealloc2d_f(par, ndim); 
    free(buf_nhat1);
    free(buf_nhat2);
    free(buf_nhat3);
    free(buf_sv1);
    free(buf_sv2);
    free(buf_sv3);
    free(buf_tp);
    free(buf_tr);
    free(buf_slip);
    free(buf_slip_kin);
    free(buf_psv);
    free(buf_slip);
    free(buf_t0);
    free(buf_tarr);
    free(buf_tarr_kin);
    free(buf_dip);
    free(buf_rake);
    free(buf_strike);
    free(full);
    free(stf);
    free(t);

    /* finalize mpi */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

float objective_function(float *p, void **extras) {
    float *svm;
    float *stf;
    float *t;
    float t0;
    int nt;
    int tmp_nt;
    int Tr;
    float slip;
    float dt;
    float res = 0.0;
    int t0_ind, tarr_ind;
    int i;

    // unpackage extras
    svm = (float*) extras[0];
    nt = *(int*) extras[1];
    dt = *(float*) extras[2];
    slip = *(float*) extras[3];
    t0_ind = *(int*) extras[4];
    tarr_ind = *(int*) extras[5];
    Tr = *(float*) extras[6];

    // if (invo_c == 0) printf("nt=%d dt=%f slip=%f t0_ind=%d tarr_ind=%d\n", nt, dt, slip, t0_ind, tarr_ind);

    t0 = t0_ind * dt;

    t = arange(0, dt*nt, dt, &tmp_nt);
    if (tmp_nt != nt) {
        printf("ERROR: sim nt and stf nt not equivalent. aborting!\n");
        exit(-1);
    }

    // compute analytical tinti function with parameter
    stf = tinti(t, nt, Tr, p[0], t0, slip);

    // compute residual as sum of squares  [ t0_ind, tarr_ind )
    for (i=t0_ind; i<tarr_ind; i++) {
        res = res + powf((svm[i] - stf[i]), 2);
    }
    free(t);
    free(stf);

    return res;
}