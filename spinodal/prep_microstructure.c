#include "spinodal.h"

float** prep_microstructure(int nx, int ny, float **conc, double c0,
double noise, float **random_ZeroToOne_array ) {
  int i, j;
//  for (i=0; i<=nx; i++){
//    for (j=0; j<=ny; j++){
//      conc[i][j] = c0 + noise*(0.5-random_ZeroToOne_array[i][j]);
//    } // end for Nx
//  } // end for Ny

    for (i = 0; i <=nx; ++i) {
        conc[i][0]=c0 + noise * (0.5 - random_ZeroToOne_array[i][0]);
    }

    for (i = 0; i <=nx; i++) {
        for (j = 0; j<=ny; ++j) {
            conc[i][j]=conc[i][0];
        }
    }
  return conc;
}
