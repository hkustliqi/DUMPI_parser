CC        = gcc
CFLAGS    = -g -O3
FFTW_ROOT = 
CPPFLAGS  = -I. -DENABLE_WARMUP 
MPICC     = mpicc

ifdef FFT_LDFLAGS
LDFLAGS  += ${FFT_LDFLAGS}
endif

ifdef COMM_LDFLAGS
LDFLAGS  += ${COMM_LDFLAGS}
endif

LIBS      = -L/usr/lib/x86_64-linux-gnu -lfftw3f -lfftw3f_threads -lm -L/home/qi/sst-dumpi/lib/ -ldumpi -lundumpi -lmpi  #-lshmem

# If PE_ENV is defined I'm assuming we're on a Cray.
ifdef PE_ENV
CC        = cc
CFLAGS    = -O3
LDFLAGS   = 
CPPFLAGS  = -I. -DENABLE_WARMUP -DNEED_SHMEM_FINALIZE
LIBS      = -lfftw3_threads -lfftw3f_threads
endif

sfft : sfft.c
	$(MPICC) $(CFLAGS) -D_MPI $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

clean :
	rm -f sfft *~

