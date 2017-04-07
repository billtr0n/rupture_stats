DBG = -Wall
FC = ftn
CC = cc
FLAGS = -O3

all:
	make mpio
	make spline
	make xapiir
	make utils
	make stats

mpio:
	$(CC) -c $(DBG) $(FLAGS) sord_mpio.c

spline:
	$(CC) -c $(DBG) $(FLAGS) spline.c

xapiir:
	$(FC) -c $(FLAGS) xapiir.f

utils:
	$(CC) -c $(DBG) $(FLAGS) utils.c
    
stats:
	$(CC) $(DBG) $(FLAGS) -o stats_mpi stats.c sord_mpio.o spline.o xapiir.o utils.o -lm -lgfortran    

clean:
	rm *.o mom_mpi

