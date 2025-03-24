CC = gcc
CFLAGS = -O3 -ffast-math
#CFLAGS = -ggdb
LIBS = -lm

all: mckaol 

mckaol: 
	$(CC) $(CFLAGS) -o mckaol mckaol.c actions.c bfsearch.c envrn.c evtlist.c futil.c lattice.c myerr.c output.c ran2.c reactions.c rxnlist.c sim.c ucell.c $(LIBS)

clean:
	rm -f *.o *~ *.Addr *.Counts *.pixie mon.out logfile
