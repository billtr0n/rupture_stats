DBG = -Wall
FC = ftn
CC = cc
FLAGS = -O3

all:
	make mpio
	make spline
	make xapiir
	make utils
	make momrate

mpio:
	$(CC) -c $(DBG) $(FLAGS) sord_mpio.c

spline:
	$(CC) -c $(DBG) $(FLAGS) spline.c

xapiir:
	$(FC) -c $(FLAGS) xapiir.f

utils:
	$(CC) -c $(DBG) $(FLAGS) utils.c
    
momrate:
	$(CC) $(DBG) $(FLAGS) -o mom_mpi ts_v1_Feb_09.c sord_mpio.o spline.o xapiir.o utils.o -lm -lgfortran    

clean:
	rm *.o mom_mpi

