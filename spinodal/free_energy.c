#include "spinodal.h"

float** free_energy(int nx, int ny, float **conc, int i, int j){
  float A = 1.0;
  // dfdcon[i][j] = A*(2.0*conc[i][j]*(1-conc[i][j])*(1-conc[i][j])
  // -2.0*(conc[i][j]*conc[i][j])*(1.0-conc[i][j]));
  // dfdcon[i][j] = A*(conc[i][j] * conc[i][j])*(1-conc[i][j])*(1-conc[i][j]);
  dfdcon[i][j] = 2.0*A*conc[i][j]*(1.0 - conc[i][j])*(1-2.0*conc[i][j]);
  return dfdcon;

}
