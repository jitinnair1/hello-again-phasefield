#include "spinodal.h"

void solve(int nx, int ny, float grad_coef,
  float dfdcon[][], float lap_con[][], lap2_con[][]){
    for(i=0; i <= nx; i++){
      for(j=0; j <= ny; j++){
        lap2_con[i][j] = dfdcon[i][j] - grad_coef*lap_con[i][j];
      }
    }
  }
