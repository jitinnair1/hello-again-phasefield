#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gperftools/profiler.h>


float** array_allocate(int nx, int ny, float **array);

void array_deallocate(int ny, float **array);

void prep_microstructure(int nx, int ny, fftw_complex *conc, double c0,
        double noise, float **random_ZeroToOne_array );

float** rand_ZeroToOne(int nx, int ny, float **random_ZeroToOne_array);

void write_to_VTK( int nx, int ny, int nz,
        double dx, double dy, double dz,
        int iprint, int NxNy, double conc[NxNy] );
