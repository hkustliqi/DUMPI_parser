NAME           = netstress
VERSION        = \"2.4.0\"
TAR_VERSION    = 2.4.0
DEFS          := -DNETSTRESS_VERSION=$(VERSION)
DEBUG          =
CC             = cc
UPC            = cc
ARCH           = generic_x86_64
SRC            = src
BUILD          = build
ifdef PE_ENV
ARCH           = CRAY
endif

ifeq "$(D)" "1"
DEBUG = -DDEBUG
endif
EXE           :=

# Compiles netstress2
all: shmem upc mpi mpi2

mpi:
	${MPICC} $(CFLAGS) -D_MPI $(DEFS) $(DEBUG) -o ./$(NAME)-mpi $(SRC)/*.c  $(LDFLAGS)

mpi2:
	${MPICC} $(CFLAGS) -D_MPI2 $(DEFS) $(DEBUG) -o ./$(NAME)-mpi2 $(SRC)/*.c $(LDFLAGS)

shmem:	
	${CC} $(CFLAGS) -D_SHMEM $(DEFS) $(DEBUG) -o ./$(NAME)-shmem $(SRC)/*.c $(LDFLAGS)
upc:	
	${UPC} -D_UPC $(DEFS) $(DEBUG) -o ./$(NAME)-upc $(SRC)/*.c $(UPCFLAGS)

clean:
	rm -f *~ *.o netstress-shmem netstress-mpi netstress-mpi2 netstress-upc

#-------------------------- Generic---------------------------------------
ifeq "$(strip $(ARCH))" "generic_x86_64"
CC          = cc
UPC         = upcc
MPICC	    = mpicc
CFLAGS      = -O3
LDFLAGS    := $(COMM_LDFLAGS) -lm -L/home/qi/sst-dumpi/lib -ldumpi -lundumpi # -lshmem
UPCFLAGS   := $(UPC_LDFLAGS) -lm
LIBS        = 
endif

#-------------------------- Cray specific modifications -----------------------
ifeq "$(strip $(ARCH))" "CRAY"
# Requires Cray programming environment
CC          = cc
UPC         = cc
MPICC	    = cc
CFLAGS      = -DCRAY -O3
UPCFLAGS    = -h upc -h omp -DCRAY -O3 -lm
LDFLAGS    := -lm
LIBS        = 
endif

#--------------------------- IBM ------------------------------------------
ifeq "$(strip $(ARCH))" "ibm_ppc64"                                             
CC          = cc                                                                
UPC         = xlupc  
MPICC	    = cc                                                           
CFLAGS      = -O3
UPCFLAGS    = -qarch=pwr7 -qtune=pwr7 -O3 -lm
LDFLAGS     = -lshmem -lm
endif        

ifeq "$(strip $(ARCH))" "ibm_ppc"
CC         := xlc
MPICC      := xlc
CFLAGS     := -qarch=pwr7 -qtune=pwr7 -q64 -O3
UPCFLAGS   := -qarch=pwr7 -qtune=pwr7 -q64 -O3 -lm
LDFLAGS     = -lm
endif 

#--------------------------- Itanium ------------------------------------------
ifeq "$(strip $ARCH))" "linux_ia64"
CC         := icc
MPICC      := icc
CFLAGS     += -ip -w1 -g
LDFLAGS    := -lm
endif 

#-------------------------- tar ----------------------------------------------
DIST_FILES = Makefile README $(SRC) $(NAME)_plot
DIST_DIRS  =
TAR_DIR    = $(NAME)-$(TAR_VERSION)
am__tar    = tar chof - "$$tardir"
tar__excl  = --exclude='*.svn' --exclude='$(SRC)/*.svn' --exclude='$(NAME)-mpi' --exclude='$(NAME)-shmem' --exclude='$(NAME)-upc'

am__remove_tardir =                                                  \
  { test ! -d $(TAR_DIR)                                             \
    || { find $(TAR_DIR) -type d ! -perm -200 -exec chmod u+w {} ';' \
         && rm -fr $(TAR_DIR); }; }


dist: 
	test -d $(TAR_DIR) || mkdir $(TAR_DIR)
	list='$(DIST_FILES)';                        \
           for file in $$list; do                    \
              echo $$file; cp -r $$file $(TAR_DIR);     \
           done;				     
	list='$(DIST_DIRS)';                         \
            for dir in $$list; do                    \
               echo $$dir; cp -a $$dir $(TAR_DIR);   \
            done

dist-gzip: dist
	tardir=$(TAR_DIR) && $(am__tar) $(tar__excl) | gzip -c >$(TAR_DIR).tar.gz
	$(am__remove_tardir)

dist-bzip2: dist
	tardir=$(TAR_DIR) && $(am__tar) $(tar__excl) | bzip2 -9 -c >$(TAR_DIR).tar.bz2
	$(am__remove_tardir)
