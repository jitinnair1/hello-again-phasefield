#include "spinodal.h"

float** solve(int nx, int ny, double grad_coef,
  float** dfdcon, float **lap_con, float **lap_dummy, int i, int j){
        lap_dummy[i][j] = dfdcon[i][j] - grad_coef*lap_con[i][j];
    return lap_dummy;
  }
