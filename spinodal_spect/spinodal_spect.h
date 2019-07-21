#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>


float** array_allocate(int nx, int ny, float **array);

fftw_complex** prep_microstructure(int nx, int ny, float **conc, float c0,
  float noise, float **random_ZeroToOne_array );

float** rand_ZeroToOne(int nx, int ny, float **random_ZeroToOne_array);

void write_to_VTK( int nx, int ny, int nz,
  float dx, float dy, float dz,
  int iprint, float **conc );
