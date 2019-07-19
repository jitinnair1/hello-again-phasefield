#include "spinodal.h"

// Generate discrete laplacian
float** laplacian(
  int nx, int ny, int nz,
  float dx, float dy, float dz,
  float **conc, float **lap_con){

    // use standard 5 point stencil
    int jp, jm, ip, im;
    int i, j;
    int count1=0, count2=0, count3=0, count4=0;

    for (i=0; i<=nx; i++ ){
      for (j=0; j<=ny; j++ ){
        jp=j+1;
        jm=j-1;

        ip=i+1;
        im=i-1;

        if(im < 0){
          im=nx;
          count1+=1;
        }

        if(ip > nx){
          ip=0;
          count2+=1;
        }

        if(jm < 0){
          jm=ny;
          count3+=1;
        }

        if(jp > ny){
          jp=0;
          count4+=1;
        }

        lap_con[i][j] =(conc[im][j] + conc[ip][j] + conc[i][jm] + conc[i][jp] -4.0*conc[i][j])/(dx*dy);
      }
    }

    // // Test of Periodic Boundary Implementation
    // printf("%d\n", count1 );
    // printf("%d\n", count2 );
    // printf("%d\n", count3 );
    // printf("%d\n", count4 );
    // printf("====\n");
    return lap_con;

  }
