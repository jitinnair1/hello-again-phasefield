#include "spinodal.h"

// Generate discrete laplacian
void laplacian(
  int nx, int ny, int nz,
  float dx, float dy, float dz,
  float conc[][], int lap_con[][]){

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
          im=Nx;
        }

        if(ip > Nx){
          ip=1;
        }

        if(jm < 0){
          jm=Ny;
        }

        if(jp > Ny){
          jp=1;
        }

        hne=con[ip][j];
        hnw=con[im][j];
        hns=con[i][jm];
        hnn=con[i][jp];
        hnc=con[i][j];

        lap_con[i][j] =(hnw + hne + hns + hnn -4.0*hnc)/(dx*dy);
      }
    }
  }
