# Compiler options
FC = gfortran
FCFLAGS = -O3 -ffast-math -funroll-loops -march=native
DEBUG = -O0 -ffast-math -funroll-loops -march=native -fcheck=bounds -g -fbacktrace
DEBUGALL = -O0 -ffast-math -funroll-loops -march=native -Wall -pedantic -Wconversion-extra -fcheck=all -g -fbacktrace

# Required libraries: Lapack, FFTW3, HDF5
LIBS = -lm -L/usr/local/lib -L/usr/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lfftw3 -llapack -lhdf5_fortran -lhdf5 -lhdf5_hl -lblas

# Source tree definitions
SRC = src
VPATH = $(SRC)
BINDIR = bin
EXEC = mueller_calc

.SUFFIXES:
.SUFFIXES: .o .mod .f90 

# Hello world
yell = "Starting Make..."

# Includes and flag for putting .mod files to directory bin, home version
INCS = -I/usr/include -I/usr/local/include/ -I/usr/include/hdf5/serial/ -J${BINDIR}

# Dependency tree
OBJECTS = ${BINDIR}/common.o \
${BINDIR}/sfunctions.o \
${BINDIR}/h5io.o \
${BINDIR}/io.o \
${BINDIR}/gaussquad.o \
${BINDIR}/integration_points.o \
${BINDIR}/translations.o \
${BINDIR}/mie.o \
${BINDIR}/possu.o \
${BINDIR}/sparse.o \
${BINDIR}/singularity_subtraction.o \
${BINDIR}/singularity_subtraction_N.o \
${BINDIR}/integrals.o \
${BINDIR}/geometry.o \
${BINDIR}/sparse_mat.o \
${BINDIR}/precorrection.o \
${BINDIR}/projection.o \
${BINDIR}/build_G.o \
${BINDIR}/gmres_module.o \
${BINDIR}/transformation_matrices.o \
${BINDIR}/setup.o \
${BINDIR}/T_matrix.o \
${BINDIR}/mueller.o \
${BINDIR}/main.o

###############################################################################

# No need to touch below, unless bad makefileing or messages need tweaking...
.PHONY: all clean
.SECONDARY: main-build

all: pre-build main-build post-build

pre-build:
	@echo $(yell)

main-build: ${EXEC} | $(BINDIR)

post-build:
	@echo "Target $(EXEC) compiled successfully"
	
debug: FCFLAGS = $(DEBUG)
debug: all

debugall: FCFLAGS = $(DEBUGALL)
debugall: all

$(BINDIR):
	@echo "Binary files are put into the directory $(BINDIR)"
	@mkdir -p ${BINDIR}

${BINDIR}/%.o: %.f90 |$(BINDIR)
	@echo "Compiling $^"
	@${FC} ${FCFLAGS} $(INCS) -c $^ -o $@ 

${EXEC}: ${OBJECTS}
	@echo "Linking the target"
	@${FC} ${FCFLAGS} ${OBJECTS} ${LIBS} -o $@

# Clean only objects
clean:
	@echo "Deleted all .o files"
	@rm -rf $(BINDIR)/*.o
	@rm -rf *.mod
	@rm -f *~

# Full clean
veryclean: clean
	@echo "Deleted all .mod files"
	@echo "Deleted executable and directory bin"
	@rm -f $(EXEC)
	@rm -rf $(BINDIR)

