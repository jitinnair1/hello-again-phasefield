#include "spinodal.h"

void free_energy(int nx, int ny, float conc[][], float dfdcon[][]){
  float A = 1.0;
  int i, j;
  for (i=0; i<=nx; i++ ){
    for(j=0; j<=ny; j++){
      dfdcon[i][j] = A*(2.0*con[i][j]*(1-con[i][j])*(1-con[i][j])
      -2.0*(con[i][j]*con[i][j])*(1.0-con[i][j]));
    }
  }

}
