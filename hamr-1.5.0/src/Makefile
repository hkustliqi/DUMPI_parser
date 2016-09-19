###############################################################################
#
#     Makefile for the Huge Asynchronous 128-bit MSB Radix sort (HAMR) benchmark.
#     Variable descriptions for the typically modified variables:
#	CC         Compiler
#       CFLAGS     Flags for the compiler (optimizations, warnings, etc...)
#       INCLUDES   Additional include paths
#       LDFLAGS    Additional paths for the linker (-L <path>)
#       LIBS       Libraries to be added
#
#     Target descriptions:
#     all:             Builds hamr
#     clean:           Removes executables and object files.
#
###############################################################################

NAME           = hamr
VERSION	       = \"1.0.0\"
TAR_VERSION    = 1.0.0

PLATFORM       = GENERIC
ifdef PE_ENV
PLATFORM       = CRAY
endif

#-------------------------- Directories ---------------------------------------
BUILD	       = build
SRC           := .
INCL          := .

SRCS          := $(notdir $(wildcard $(SRC)/*.c))

#OBJS           = $(patsubst %.c,$(BUILD)/%.o,$(notdir $(filter %.c,$(SRC))))
OBJS           = $(patsubst %.c,$(BUILD)/%.o,$(SRCS))

#-------------------------- Executables ---------------------------------------
EXE	      := $(NAME)

# Dependency compiler
MAKEDEPC    = gcc -MM

CC          = mpicc #gcc
CFLAGS      = -O3 -std=c99
INCFLAGS   := $(patsubst %,-I%,$(INCL))
DEFS       := -DPACKAGE_VERSION=$(VERSION) -D$(PLATFORM) -D_SGI_SOURCE
LDFLAGS    :=
LIBS        = -lm -L/home/qi/sst-dumpi/lib -ldumpi -lundumpi -lmpi #-lsma -lm 

#-------------------------- Cray specific modifications -----------------------
ifeq "$(strip $(PLATFORM))" "CRAY"
# Requires Cray programming environment
CC          = cc
CFLAGS      = -O3 -g
DEFS       := -DPACKAGE_VERSION=$(VERSION) -D$(PLATFORM)
LDFLAGS    := 
LIBS        = 
ifeq "$(strip $(PE_ENV))" "GNU"
CFLAGS     += -std=c99 
endif
endif

#--------------------------- Misc ---------------------------------------------
.PHONY: setup dist clean cleanall depclean doc

.SUFFIXES:
.SUFFIXES: .o .c .h .a

# places to search for various files
vpath %.c $(SRC)
vpath %.h $(INCL)

all:   $(EXE) 

#--------------------------- Executables -------------------------------------
$(EXE): $(OBJS)
	$(CC) -o $@ $(LDFLAGS) $^ $(LIBS) 

#-------------------------- Object file creation -----------------------------
ifeq ($(strip $(filter dist clean%, $(MAKECMDGOALS)) ),)
-include $(OBJS:.o=.d)
endif

$(BUILD)/%.d: %.c
	@if [ ! -d $(BUILD)  ]; then mkdir -p $(BUILD); fi 
	@echo "Generating dependencies for $<"
	@$(MAKEDEPC) $(DEFS) $(INCFLAGS) $< | \
        sed 's@^\(.*\)\.o:@$(BUILD)/\1.d $(BUILD)/\1.o:@g' > $@

$(BUILD)/%.o: %.c
	@if [ ! -d $(BUILD)  ]; then mkdir -p $(BUILD); fi
	$(CC) -o $@ -c $(CFLAGS) $(INCFLAGS) $(DEFS) $< 

#-------------------------- tar ----------------------------------------------
DIST_FILES = Makefile README *.c *.h
DIST_DIRS  = 
TAR_DIR    = $(NAME)_$(TAR_VERSION)
am__tar    = tar chof - "$$tardir"
tar__excl  = --exclude='*.svn' --exclude $(BUILD) --exclude $(EXE)

am__remove_tardir =                                                  \
  { test ! -d $(TAR_DIR)                                             \
    || { find $(TAR_DIR) -type d ! -perm -200 -exec chmod u+w {} ';' \
         && rm -fr $(TAR_DIR); }; }


dist: 
	test -d $(TAR_DIR) || mkdir $(TAR_DIR)
	list='$(DIST_FILES)';                        \
           for file in $$list; do                    \
              echo $$file; cp $$file $(TAR_DIR);     \
           done
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

#-------------------------- documentation ------------------------------------
doc:
	doxygen Doxyfile
	cd doc/latex && make && pdflatex refman


#-------------------------- clean up -----------------------------------------
cleanall: clean depclean

clean:
	@echo Cleaning entire distribution...
	@rm -rf $(BUILD)/*.o $(EXE) $(LIB)

depclean: 
	@rm -rf $(BUILD)/*.d
