#include "spinodal.h"

float** solve(int nx, int ny, float grad_coef,
  float** dfdcon, float **lap_con, float **lap_dummy){
    int i, j;
    for(i=0; i < nx; i++){
      for(j=0; j < ny; j++){
        lap_dummy[i][j] = dfdcon[i][j] - grad_coef*lap_con[i][j];
      }
    }
    return lap_dummy;
  }
