
######################
## profstrchim_gi (general input with ASAs) Makefile ##
######################

## Comment the following to hide warnings
WARNINGS="-Wall -Wextra"

OPTS="-O2  -pedantic"
#CODS=ConfigFile.cpp
EXE=WSME_genTden_loopy

## libraries
LIB="-lm"

gfortran $WARNINGS $OPTS -c parameter_modules.f90

gfortran $WARNINGS $OPTS -c calc_e_Phi_genTden_loopy.f90

gfortran $WARNINGS $OPTS -c calc_thermoab.f90

gfortran $WARNINGS $OPTS -c calc_thermo_loopy.f90

gfortran $WARNINGS $OPTS -c stringhe.f90

gfortran $WARNINGS $OPTS -c disulfide_module.f90

gfortran $WARNINGS $OPTS -c WSME_genTden_loopy.f90

gfortran -o $EXE parameter_modules.o calc_e_Phi_genTden_loopy.o calc_thermoab.o calc_thermo_loopy.o  stringhe.o disulfide_module.o WSME_genTden_loopy.o ${LIB};

mv *.mod Modules/

rm *.o

rm *.ent