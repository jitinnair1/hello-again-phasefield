#include "spinodal.h"

float** solve2(int nx, int ny, double dt,
  double mobility, float **lap2_con, float **conc, int i, int j){
    conc[i][j] = conc[i][j] + dt*mobility*lap2_con[i][j];
    return conc;
  }
