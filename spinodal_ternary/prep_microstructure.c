#include "ternary_coarsening.h"

void prep_microstructure(int nx, int ny, fftw_complex *conc, float c0,
float noise, float **random_ZeroToOne_array ) {
  int ii=0;
  int i, j;
  for (i=0; i < nx; i++){
    for (j=0; j < ny; j++){
      ii=i*nx+j;
      conc[ii] = c0 + noise*(0.5-random_ZeroToOne_array[i][j]);
    } // end for Nx
  } // end for Ny
}
