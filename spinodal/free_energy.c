#include "spinodal.h"

float** free_energy(int nx, int ny, float **conc, float **dfdcon){
  float A = 1.0;
  int i, j;
  for (i=0; i<=nx; i++ ){
    for(j=0; j<=ny; j++){
      // dfdcon[i][j] = A*(2.0*conc[i][j]*(1-conc[i][j])*(1-conc[i][j])
      // -2.0*(conc[i][j]*conc[i][j])*(1.0-conc[i][j]));
      dfdcon[i][j] = A*(conc[i][j] * conc[i][j])*(1-conc[i][j])*(1-conc[i][j]);
    }
  }
  return dfdcon;

}
