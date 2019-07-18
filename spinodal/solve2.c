#include "spinodal.h"

void solve2(int nx, int ny, float dt,
  float mobility, float lap2_con[nx][ny], float **conc){
    int i, j;
    for (i=0; i<nx; i++){
      for (j=0; j<ny; j++){
        conc[i][j] = conc[i][j] + dt*mobility*lap2_con[i][j];
      } // end for Nx
    } // end for Ny

}
