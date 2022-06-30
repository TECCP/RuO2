###############################################################################
#                                                                             #
#                          IBM RS6000  AIX 3.2.x                              #
#                                                                             #
#  make           ... generate executable for the REAL-version                #
#  make complex   ... generate executable for the COMPLEX-version             #
#  make clean     ... delete unnecessary files                                #
#                                                                             #
###############################################################################
#
# FC ........... compiler name
# FOPT ......... optimization flags
# FGEN ......... code generation flags (flags not related to optimization)
# LDFLAGS ...... linker flags
# R_LIBS ....... libraries needed to build the REAL executable
# C_LIBS ....... libraries needed to build the COMPLEX executable
# DESTDIR ...... dir. where the executable should go (without trailing '/'!)
# R_EXECNAME ... name of the resulting REAL executable (without prefixed path!)
# C_EXECNAME ... name of the resulting COMPLEX executable ( -"- )
#
SHELL = /bin/sh
FC = mpif77
CC = icc
FOPT = -O2 -funroll-loops -assume byterecl -align -pad -static-intel
FGEN = 
LDFLAGS = 
R_LIBS = -L/opt/intel-14.0.3/composer_xe_2013_sp1.3.174/mkl/lib/intel64 \
-lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread
DESTDIR	= .
R_EXECNAME = dmft

###############################################################################

TYPE='REAL'
TYPE_COMMENT='!_REAL'
FFLAGS = $(FGEN) $(FOPT)
R_EXEC = $(DESTDIR)/$(R_EXECNAME)
C_EXEC = $(DESTDIR)/$(C_EXECNAME)

#..............................................................................
#
#  Object files common to both REAL and COMPLEX type subroutines
#
#
OBJS    = main.o xn.o mu_calc.o  gloc.o w2t1.o energy.o asym.o spec.o cputim.o \
          delta_as.f nmat.o delta.o mu_calc_newton.o xnd.o


#..............................................................................
#
#  Object files for REAL and COMPLEX versions
#
R_OBJS = $(OBJS) 
C_OBJS = $(OBJS)

#..............................................................................
#
#  Some source files contain both REAL and COMPLEX versions. These files
#  have the file extension '.frc' and must be treated differently.
#
.SUFFIXES:	.frc

#..............................................................................
#
#  Build executable (either REAL or COMPLEX versions)
#
real: rm_c_files
	$(MAKE) TYPE='REAL' TYPE_COMMENT='!_REAL' $(R_EXEC)

complex: rm_r_files
	$(MAKE) TYPE='COMPLEX' TYPE_COMMENT='!_COMPLEX' $(C_EXEC)

#..............................................................................

$(R_EXEC): $(R_OBJS)
	$(FC) -o $(R_EXEC) $(R_OBJS) $(LDFLAGS) $(R_LIBS)
	size $(R_EXEC)

$(C_EXEC): $(C_OBJS)
	$(FC) -o $(C_EXEC) $(C_OBJS) $(LDFLAGS) $(C_LIBS)
	size $(C_EXEC)

rm_c_files:
	if [ -f .complex ]; then \
	   $(MAKE) clean; \
	fi
	touch .real

rm_r_files:
	if [ -f .real ]; then \
	   $(MAKE) clean; \
	fi
	touch .complex

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#
$(R_OBJS):		
$(C_OBJS):		

#..............................................................................
#
#  remove unnecessary files (executable(s) are not removed)
#
# Dependancies
$(OBJS):		Makefile

clean:
	-rm -f *.o _tmp_.* *.P .real .complex *.mod

#..............................................................................
#
#  define inference rules to generate object files from source files
#  (From .frc files the appropriatly typed subroutine (REAL/COMPLEX) has to
#  be extracted first.)
#
.f.o:
	$(FC) $(FFLAGS) -c $<
.frc.o:
	@sed -e 's/^$(TYPE_COMMENT)//'   $< >$*.f
	@echo "$<: $(TYPE) version extracted"
	$(FC) $(FFLAGS) -c $*.f
	rm $*.f
