# Use bash as the shell for evaluating expressions
SHELL=/bin/bash

# Define the source code files
SRC=$(shell echo src/*.c)

# Define header files (mostly for the backup)
HDR=$(shell echo src/*.h)

# Shell substitution turns source file names into object file names
OBJ=$(subst .c,.o,${SRC})
#########OBJECTS = $(SOURCES:.c=.o)#another way from internet

# Get the current version of the compiler from
# /usr/bin/gcc -v 2>&1 | grep -w version | cut -d ' ' -f3
# and use that number here
GCCV=5.4.0

ARCH=x86_64-linux-gnu
# ARCH=x86_64-redhat-linux on Fedora machines (usually)

# Include paths that should, but don't always, appear automagically
INCLUDEPATHS=-I/usr/include/c++/${GCCV} \
	-I/usr/include/c++/${GCCV}/${ARCH} \
	-I/usr/lib/gcc/${ARCH}/${GCCV}/include \
	-I/usr/include/gsl \
        -I/usr/include/${ARCH} \
	-I/usr/include/x86_64-linux-gnu/c++/5 \
	-I/home/mcdermott/Codes-Scripts/Vortex_2D/src  #find the header dammit


# Source files to be compiled as a library
#LIBSRC=parameter.cpp

# Substitution for the object file names
#LIBOBJ=$(#subst .cpp,.o,${LIBSRC})

# Name of the output library
#LIBNAME=parameter

# Compiler directive for including the library in the program
#LIBRARIES=-l${LIBNAME}

# Where to look for the library object
#LIBSEARCH=-L./

# Compiler options
DEBUG=0
ifeq (${DEBUG},1)
  CPUOPT=-g3 -Wall -Wunused -pg -fno-strict-aliasing -finline-functions -DDEBUG
  LIBS=${LIBRARIES} -pg -lm -lgsl -lgslcblas -std=c++11 -L/usr/lib/x86_64-linux-gnu -L/usr/include/x86_64-linux-gnu/c++/5/
else
  CPUOPT= -O3 -mhard-float -msse2 -fno-strict-aliasing
  LIBS=${LIBRARIES} -pg -lm -lgsl -lgslcblas #-L/usr/lib/x86_64-linux-gnu 
endif

# compile command
CC=gcc $(CPUOPT) $(INCLUDEPATHS)

# Link command
LINK=gcc -o $(BIN) $(LIBS) #$(LIBSEARCH) $(LIBS)

# Name of the executable 
BIN=vortex_geometry

# Make sections... call e.g. "make clean" to run the commands after clean:

### compile and link everything
#all:	dep lib $(BIN) #commented out because of "lib"
all:	dep $(BIN)
# Clean up & remove unwanted files
clean:
	rm -f $(BIN) $(OBJ) $(LIBOBJ) $(LIBNAME) *~ *.bak .*.bak gmon.out

# Clean up but leave backup files
tidy:
	rm -f $(BIN) $(OBJ) $(LIBOBJ) $(LIBNAME) 

# Force everything to recompile
force:	tidy all

# If any object file has changed, re-link the binary object
${BIN}:	dep $(OBJ)
	@echo -e "\n>>>>>>>>>>>> Linking <<<<<<<<<<<<<"
	$(LINK) $(OBJ) $(LIBS) #double the libs!
        ifeq (${DEBUG},0)
	  @strip ${BIN}
        endif

# If any source file has changed, recompile that source file
.c.o:
	@echo -e "\n>>>>>>>>>>>> Compiling $< -> $@ <<<<<<<<<<<<<"
	$(CC) -c $< -o $@

### If the library file(s) have changed, remake the library
#lib:
#	@echo -e "\n>>>>>>>>>>>> Making ${LIBNAME} Library <<<<<<<<<<<<<"
#	$(CC) -fPIC $(CPUOPT) -c ${LIBSRC}
#	$(CC) -shared -o lib${LIBNAME}.so ${LIBOBJ}

# Make a tarball of important files
backup:
	@tar -zcf ${BIN}-backup.tar.gz $(SRC) $(LIBSRC) $(HDR) \
	.dependencies Makefile

# Get the system dependencies and save them
depend dep:
	makedepend $(INCLUDEPATHS) $(SRC) -f .dependencies;

# Include all the system dependencies in case any of them change.
# Caveat: Make reads in this file, and tries to include .dependencies
# *before* doing anything else. So there's no way to create it here
# if it doesn't exist before make tries to include it. Stoopid make.
# So issue a "touch .dependencies" from the command line if it doesn't
# exists before doing a "make", otherwise make will choke without 
# being able to do a makedepend to create the file.
include .dependencies 
