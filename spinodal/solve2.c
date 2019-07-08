#include "spinodal.h"

void solve2(int nx, int ny, float dt,
  float mobility, float lap2_con, float conc){

    for (i=0; i<=nx; i++){
      for (j=0; j<=ny; j++){
        conc[i][j] = con[i][j] + dtime*mobility*lap2_con[i][j];
      } // end for Nx
    } // end for Ny

}
