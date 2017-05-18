FC = mpixlf90
CC = mpixlc
FLAGS = -g -O3

all:
	make sord_mpio
	make utils
	make stf
	make sord_utils
	make brute
	make stats

sord_mpio:
	$(CC) -c $(FLAGS) sord_mpio.c

utils:
	$(CC) -c $(FLAGS) utils.c
    
sord_utils:
	$(CC) -c $(FLAGS) sord_utils.c

stf:
	$(CC) -c $(FLAGS) stf.c

brute:
	$(CC) -c $(FLAGS) brute.c

stats:
	$(CC) $(FLAGS) -o stats_mpi stats.c sord_mpio.o utils.o sord_utils.o brute.o stf.o -lm
	

clean:
	rm *.o stats_mpi

