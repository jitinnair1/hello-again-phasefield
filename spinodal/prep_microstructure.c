#include "spinodal.h"

void prep_microstructure(int nx, int ny, float **conc, float c0,
float noise, float **random_ZeroToOne_array ) {
  int i, j;
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      conc[i][j] = c0 + noise*(0.5-random_ZeroToOne_array[i][j]);
    } // end for Nx
  } // end for Ny
}
