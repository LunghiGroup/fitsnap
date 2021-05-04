export LAMMPS_DIR=/home/Alessandro/Documents/Sources/lammps_git/lammps

mpif90 -O2 -g -c  *.f90  -I${LAMMPS_DIR}/examples/COUPLE/fortran2 

mpif90 *.o -L${LAMMPS_DIR}/src -L${LAMMPS_DIR}/examples/COUPLE/fortran2 -llammps_fortran -llammps_mpi -llapack -lmpi_cxx -lstdc++ -lm -o fitsnap.x
mv fitsnap.x bin/
