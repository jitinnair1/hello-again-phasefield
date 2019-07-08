#include "spinodal.h"

// Generate discrete laplacian
void laplacian(
  int nx, int ny, int nz,
  float dx, float dy, float dz,
  float conc[nx][ny], float lap_con[nx][ny]){

    // use standard 5 point stencil
    int jp, jm, ip, im;
    int i, j;

    for (i=0; i<nx; i++ ){
      for (j=0; j<ny; j++ ){
        jp=j+1;
        jm=j-1;

        ip=i+1;
        im=i-1;

        if(im < 0){
          im=nx;
        }

        if(ip > nx){
          ip=1;
        }

        if(jm < 0){
          jm=ny;
        }

        if(jp > ny){
          jp=1;
        }

        lap_con[i][j] =(conc[im][j] + conc[ip][j] + conc[i][jm] + conc[i][jp] -4.0*conc[i][j])/(dx*dy);
      }
    }
  }
