#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gperftools/profiler.h>

void prep_microstructure(int iflag, int Nx, int Ny, int NxNy, fftw_complex conc[NxNy], double conc_print[NxNy], double c0,
                         double noise, double random_ZeroToOne_array[NxNy] ) ;

double *rand_ZeroToOne(int Nx, int Ny, int seed, double *random_ZeroToOne_array) ;

void write_to_VTK( int nx, int ny, int nz,
        double dx, double dy, double dz,
        int iprint, int NxNy, double conc[NxNy] );

void write_init_conc(int Nx, int Ny, int NxNy, double conc_print[NxNy]);
