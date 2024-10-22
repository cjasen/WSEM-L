######################
## profstrchim_gi (general input with ASAs) Makefile ##
######################

## Comment the following to hide warnings
WARNINGS="-Wall -Wextra"

OPTS="-O2 -pedantic"
#CODS=ConfigFile.cpp
EXE=optimization_WSME_loopy

## libraries
LIB="-lm"

gfortran $WARNINGS $OPTS -c parameter_modules.f90

gfortran $WARNINGS $OPTS -c optimization_WSME_loopy.f90

gfortran -o $EXE parameter_modules.o optimization_WSME_loopy.o ${LIB} -L"C:/Users/cjase/gcc/nlopt/lib" -lnlopt

#mv *.mod Modules/

rm *.o
