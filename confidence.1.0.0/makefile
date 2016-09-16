# Confidence Makefile
#
# where are the MPI libraries and which to include?
#    MPI_LIBDIR = -L/usr/lib64/openmpi/1.3.1-gcc
#    MPI_LIBS   = -lmpi
MPI_LIBDIR = 
MPI_LIBS   = -lmpi

# compiler and flags
CC       =  mpicc #gcc
CFLAGS   = -m64 -march=core2 -mtune=core2 -g -O3 -lm -L/home/qi/sst-dumpi/lib -ldumpi -lundumpi

# libraries
LIBFLAGS = $(MPI_LIBDIR) $(MPI_LIBS)

COLLATERAL = LICENSING.pdf ACKNOWLEDGEMENTS COPYRIGHT FAQ ORBTIMER scripts/config.sh makefile
SRCS     = confidence.c comm.c comm_test1.c comm_test2.c comm_test3.c measurement.c options.c orbtimer.c
HDRS     = confidence.h orbtimer.h
OBJS     = confidence.o comm.o comm_test1.o comm_test2.o comm_test3.o measurement.o options.o orbtimer.o

confidence    : $(OBJS)
	$(CC)  -o confidence $(LIBFLAGS) $(OBJS) $(CFLAGS)
	cp confidence ./bin/confidence

orbtimer.o    : config.h orbtimer.h orbtimer.c

confidence.o  : config.h confidence.h confidence.c

measurement.o : config.h confidence.h measurement.c

options.o     : config.h confidence.h options.c

comm.o        : config.h confidence.h comm.c

comm_test1.o  : config.h confidence.h comm_test1.c

comm_test2.o  : config.h confidence.h comm_test2.c

comm_test3.o  : config.h confidence.h comm_test3.c

config.h      :
	echo No config.h
	echo you need to run scripts/config.sh 
	scripts/config.sh

LICENSING.pdf : LICENSING.doc
	oowriter LICENSING.doc

dist          : $(COLLATERAL) $(SRCS) $(HDRS)
	scripts/build_release

commit        :
	git commit -a

clean         : 
	rm -f *.o confidence 

cleanall       : clean 
	rm -f ./bin/confidence

distclean     :
	rm -f *.o confidence config.h


edit          :
	gvim -p *.h makefile TODO
	gvim -p confidence.c measurement.c options.c orbtimer.c
	gvim -p comm*.c

run:
	mpirun 2 ./bin/confidence
