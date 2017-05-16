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
#include "stf.h"

#define ndim 1
#define DEFAULT -99999.0f

// function prototypes
float objective_function(float *p, void **extras);

int main (int argc, char*argv[]) {
    //////////////////////////////////////////////////////////////////////////
    // USER PARAMETERS
    //////////////////////////////////////////////////////////////////////////
    int nchunks = 1; /* also equal to the number of i/o calls */
    int nx = 2601, ny = 801, nt=10001;
    float dt = 0.002;
    int dx = 25;
    float pl[ndim], ph[ndim], pd[ndim];
    int debug = 0;
    int nfields = 12;
    
    // input file names
    char *nhat1_file = "./out/nhat1";
    char *nhat2_file = "./out/nhat2";
    char *nhat3_file = "./out/nhat3";
    char *sv1_file = "./out/sv1";
    char *sv2_file = "./out/sv2";
    char *sv3_file = "./out/sv3";

    // output file names
    char *tp_file = "./stats/tp";
    char *tr_file = "./stats/tr";
    char *slip_file = "./stats/slip";
    char *slip_kin_file = "./stats/slip_kin";
    char *psv_file = "./stats/psv";
    char *psv_kin_file = "./stats/psv_kin";
    char *t0_file = "./stats/t0";
    char *tarr_file = "./stats/tarr";
    char *tarr_kin_file = "./stats/tarr_kin";
    char *dip_file = "./stats/dip";
    char *rake_file = "./stats/rake";
    char *strike_file = "./stats/strike";
    
    //////////////////////////////////////////////////////////////////////////
    // Don't modify anything below here. 
    //////////////////////////////////////////////////////////////////////////
    obj_fun_ptr func = &objective_function;

    float *buf_sv1, *buf_sv2, *buf_sv3;
    float **sv1, **sv2, **sv3, **svm;
    float *buf_nhat1, *buf_nhat2, *buf_nhat3;

    float *buf_tp, *buf_tr;
    float *buf_tarr_kin;
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
    float *min;
    float *full;

    float *rbuf;

    int tarr_ind, tarr_kin_ind;
    int tmp_n;
    int n_comb;
    int *n_par;
    int t0_ind;

    void *extras[25];

    int csize;
    int rank, nprocs;
    int s0, xi, yi;
    int k, n, l, i;

    double t1, t2;
    
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
    buf_sv1 = calloc(csize*nt, sizeof(float));
    buf_sv2 = calloc(csize*nt, sizeof(float));
    buf_sv3 = calloc(csize*nt, sizeof(float));
    buf_nhat1 = calloc(csize, sizeof(float));    
    buf_nhat2 = calloc(csize, sizeof(float));
    buf_nhat3 = calloc(csize, sizeof(float));

    /* 2D arrays for time-series */
    sv1 = (float**)calloc(csize, sizeof(float*));
    sv2 = (float**)calloc(csize, sizeof(float*));
    sv3 = (float**)calloc(csize, sizeof(float*));
    svm = (float**)calloc(csize, sizeof(float*));
    
    for (l=0; l<csize; l++) {
        sv1[l] = (float*) calloc(nt, sizeof(float));
        sv2[l] = (float*) calloc(nt, sizeof(float));
        sv3[l] = (float*) calloc(nt, sizeof(float));
        svm[l] = (float*) calloc(nt, sizeof(float));
    }

    // calculated arrays
    buf_tp = calloc(csize, sizeof(float));
    buf_tr = calloc(csize, sizeof(float));
    buf_strike = calloc(csize, sizeof(float));
    buf_dip = calloc(csize, sizeof(float));
    buf_rake = calloc(csize, sizeof(float));
    buf_slip = calloc(csize, sizeof(float));
    buf_slip_kin = calloc(csize, sizeof(float));
    buf_psv = calloc(csize, sizeof(float));
    buf_psv_kin = calloc(csize, sizeof(float));
    buf_t0 = calloc(csize, sizeof(float));
    buf_tarr = calloc(csize, sizeof(float));
    buf_tarr_kin = calloc(csize, sizeof(float));

    n_par = malloc(sizeof(int) * ndim);
    par = malloc(sizeof(float*) * ndim);


    // allocate buffers for receiving data from all processes
    if (rank==0) {
        rbuf = malloc(sizeof(float)*csize*nprocs);
    }

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
            fprintf(stderr, "rank=%i s0=%d off=%lld xi=%d yi=%d\n", rank, s0, off, xi, yi);
        }
                
        /* read dynamic rupture parameters */ 
        if (rank == 0) fprintf(stdout, "Reading time-series...\n");
        MPI_Barrier(MPI_COMM_WORLD);
        t1 = MPI_Wtime();
        read_time_series(sv1_file, nx, ny, xi, yi, nt, csize, buf_sv1);
        read_time_series(sv2_file, nx, ny, xi, yi, nt, csize, buf_sv2);
        read_time_series(sv3_file, nx, ny, xi, yi, nt, csize, buf_sv3);
        MPI_Barrier(MPI_COMM_WORLD);
        t2 = MPI_Wtime();
        if (rank == 0) fprintf(stdout, "Finished reading in %f seconds.\n", t2-t1);

        /* these routines might not need to be mpi/io files are quite small
           might be better to read in at rank zero process and broadcast */
        if (rank == 0) fprintf(stdout, "Reading fault parameters...\n");
        MPI_Barrier(MPI_COMM_WORLD);
        t1 = MPI_Wtime();
        read_fault_params(nhat1_file, off, csize, buf_nhat1);
        read_fault_params(nhat2_file, off, csize, buf_nhat2);
        read_fault_params(nhat3_file, off, csize, buf_nhat3);
        MPI_Barrier(MPI_COMM_WORLD);
        t2 = MPI_Wtime();
        if (rank == 0) fprintf(stdout, "Finished reading in %f seconds.\n", t2-t1);

        /* organize into 2d arrays for easy filtering */
        for (n=0; n < nt; n++) {
            for (l=0; l < csize; l++) {
                sv1[l][n] = buf_sv1[n*csize+l];
                sv2[l][n] = buf_sv2[n*csize+l];
                sv3[l][n] = buf_sv3[n*csize+l];
                svm[l][n] = sqrt(powf(sv1[l][n],2)+powf(sv2[l][n],2)+powf(sv3[l][n],2));
            }
        }

        // perform evalulations of slip-rate files

        MPI_Barrier(MPI_COMM_WORLD);
        if (rank==0) fprintf(stdout, "Calculating slip-rate statistics for chunk %d.\n", k);
        t1 = MPI_Wtime();
        for (l=0; l<csize; l++) {
            if ((rank==0) && (l % 100==0)) fprintf(stdout, "Finished processing subfault %d.\n", l);
            // compute peak slip velocity magnitude
            buf_psv[l] = maximum( svm[l], nt );
            
            // compute strike, dip, rank only if patch slipped.
            buf_strike[l] = get_strike( buf_nhat1[l], buf_nhat3[l] );
            buf_dip[l] = get_dip( buf_nhat1[l], buf_nhat2[l], buf_nhat3[l] );

            // assemble dynamic slip vector and normalize for angle estimations
            su1 = trapz( sv1[l], nt, dt );
            su2 = trapz( sv2[l], nt, dt );
            su3 = trapz( sv3[l], nt, dt );
            su = unit_vector_3d( su1, su2, su3 );
            // if (rank == 4) fprintf(stdout, "(%d:%d) buf_psv[l]=%f su1=%f su2=%f su3=%f\n", rank, l, buf_psv[l], su1, su2, su3);
            
            // ignore subfaults that did not slip, and assign default values
            if (buf_psv[l] <= 0.001) {
                // fprintf(stdout, "(%d) max(svm) = %f on fault node. Assigning default values.\n", rank);
                buf_t0[l] = DEFAULT;
                buf_tarr_kin[l] = DEFAULT;
                buf_tarr[l] = DEFAULT;
                buf_tr[l] = DEFAULT;
                buf_tp[l] = DEFAULT;
                buf_psv_kin[l] = DEFAULT;
                buf_slip_kin[l] = DEFAULT;
                buf_rake[l] = DEFAULT;

                // should be quite close to zero.
                buf_slip[l] = sqrt(powf(su1,2)+powf(su2,2)+powf(su3,2));

            // subfault ruptured, svm > 0.001
            } else {
                t0_ind = find_first( svm[l], nt, ">", 0.001 );
                buf_t0[l] = t0_ind * dt;

                // compute kinematic rupture arrest
                // note: defined as first occurence after t0 where svm drops below 0.001
                tarr_kin_ind = t0_ind + find_first( &svm[l][t0_ind+1], nt-t0_ind-1, "<", 0.001 );
                buf_tarr_kin[l] = tarr_kin_ind * dt;

                // compute dynamic rupture arrest
                // note: defined where 
                tarr_ind = find_last( svm[l], nt, ">", 0.001 );
                if (tarr_ind == -1) { 
                    buf_tarr[l] = DEFAULT;
                } else {
                    buf_tarr[l] = tarr_ind * dt;
                }

                // compute rise time
                buf_tr[l] = buf_tarr_kin[l] - buf_t0[l];
                // if (rank==4) fprintf(stdout, "(%d:%d) t0_ind=%d tarr_kin_ind=%d t0=%f tarr_kin=%f Tr=%f\n", 
                                     //rank, l, t0_ind, tarr_kin_ind, buf_t0[l], buf_tarr_kin[l], buf_tr[l]);

                // compute kinematic slip
                tmp_n = tarr_kin_ind - t0_ind + 1;

                // if (rank==0) fprintf(stdout, "(%d:%d) tmp_n=%d\n", rank, l, tmp_n);
                buf_slip_kin[l] = trapz( &svm[l][t0_ind], tmp_n, dt );

                // compute dynamic slip
                buf_slip[l] = trapz( svm[l], nt, dt );

                // compute peak slip velocity
                buf_psv[l] = maximum( svm[l], nt );
                // if (rank==4) fprintf(stdout, "(%d:%d) k_slip=%f d_slip=%f psv=%f\n", rank, l, buf_slip_kin[l], buf_slip[l], buf_psv[l]);
                
                // compute rake
                buf_rake[l] = get_rake( buf_nhat1[l], buf_nhat2[l], buf_nhat3[l], su[0], su[1], su[2] );

                // here 50 is arbitrary, but we want "sufficiently" long ruptures for fitting.
                // note: 50*dt = 0.1 s
                if (buf_tr[l] > 50*dt) {

                    // set-up extra information for objective function
                    extras[0] = (void*) svm[l];
                    extras[1] = (void*) &nt;
                    extras[2] = (void*) &dt;
                    extras[3] = (void*) &buf_slip_kin[l];
                    extras[4] = (void*) &t0_ind;
                    extras[5] = (void*) &tarr_kin_ind;
                    extras[6] = (void*) &buf_tr[l];

                    // if (rank==0) fprintf(stderr, "n_comb=%d ndim=%d\n", n_comb, ndim);

                    // define parameters for estimation
                    // same as: arange(dt, Tr/2.0, dt);
                    pl[0] = dt;
                    ph[0] = buf_tr[l]/2.0 - dt;
                    pd[0] = dt;

                    // create possibly jagged array of parameters       
                    for (i=0; i<ndim; i++)  {
                        par[i] = arange(pl[i], ph[i], pd[i], &n_par[i]);
                    }

                    // generate all possible combinations of parameters
                    comb = cartesian_product( par, n_par, &n_comb, ndim );

                    // find best fitting peak-time
                    min = fmin_brute( func, comb, ndim, n_comb, extras, NULL );
                    buf_tp[l] = min[0];
                    // fprintf(stderr, "(%d:%d) Tp=%f Tr=%f\n", rank, l, buf_tp[l], buf_tr[l]);

                    // compute tinti function with best fitting parameters
                    stf = tinti(t, nt, buf_tr[l], min[0], buf_t0[l], buf_slip_kin[l]);

                    // compute kinematic psv
                    buf_psv_kin[l] = maximum( stf, nt );            

                    // free memory allocated during loop
                    free(min);
                    free(comb);
                    free(stf);
                    
                // rupture insufficiently long assuming default parameters for fitted fields
                } else {
                    buf_tp[l] = DEFAULT;
                    buf_psv_kin[l] = DEFAULT;
                }
            }
            // if (rank==0) fprintf(stdout, "strike=%f dip=%f rake=%f\n", buf_strike[l], buf_dip[l], buf_rake[l]);

            // end nicely.
            free(su);

        } /*end csize loop */
        // MPI_Barrier(MPI_COMM_WORLD);
        t2 = MPI_Wtime();
        fprintf(stdout, "(%d) Finished calculating chunk in %f seconds.\n", rank, t2-t1);
            
        /* write out necessary info */
        if (rank==0) fprintf(stdout, "Preparing to gather calculated statistics.\n");
        MPI_Barrier(MPI_COMM_WORLD);
        t1 = MPI_Wtime();
        MPI_Gather(buf_tp, csize, MPI_FLOAT, rbuf, csize, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        t2 = MPI_Wtime();
        if (rank==0) fprintf(stdout, "Receiving processed data using MPI_Gather() in %f seconds.\n", t2-t1);
        
        //write_fault_params(tp_file, off, csize, buf_tp, MPI_FLOAT); 
        //write_fault_params(tr_file, off, csize, buf_tr, MPI_FLOAT); 
        //write_fault_params(slip_file, off, csize, buf_slip, MPI_FLOAT);
        //write_fault_params(slip_kin_file, off, csize, buf_slip_kin, MPI_FLOAT); 
        //write_fault_params(psv_file, off, csize, buf_psv, MPI_FLOAT); 
        //write_fault_params(psv_kin_file, off, csize, buf_psv_kin, MPI_FLOAT); 
        //write_fault_params(t0_file, off, csize, buf_t0, MPI_FLOAT); 
        //write_fault_params(tarr_file, off, csize, buf_tarr, MPI_FLOAT); 
        //write_fault_params(tarr_kin_file, off, csize, buf_tarr_kin, MPI_FLOAT);
        //write_fault_params(dip_file, off, csize, buf_dip, MPI_FLOAT);
        //write_fault_params(rake_file, off, csize, buf_rake, MPI_FLOAT);
        //write_fault_params(strike_file, off, csize, buf_strike, MPI_FLOAT);
        
        // just write out right now, later, we will need to handle the case
        // when nchunks is not equal to 1.
        if (rank==0) {
            t1 = MPI_Wtime();
            FILE *fout = fopen("./stats/tp", "wb");
            fwrite(rbuf, sizeof(float)*csize*nprocs, 1, fout);
            fclose(fout);
            t2 = MPI_Wtime();
            fprintf(stdout, "Finished writing in %f seconds.\n", t2-t1);
        }

    } /* end main loop */


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
    free(rbuf);
    free(buf_tarr);
    free(buf_tarr_kin);
    free(buf_dip);
    free(buf_rake);
    free(buf_strike);
    //free(full);
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
    float Tr;
    float slip;
    float dt;
    float res = 0.0;

    int nt;
    int tmp_nt;
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

    // fprintf(stdout, "(in obj fun): nt=%d dt=%f slip=%f t0_ind=%d tarr_ind=%d Tr=%f\n", nt, dt, slip, t0_ind, tarr_ind, Tr);

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
