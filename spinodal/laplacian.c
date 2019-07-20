#include "spinodal.h"

// Generate discrete laplacian
float** laplacian(
  int nx, int ny, int nz,
  float dx, float dy, float dz,
  float **conc, int i, int j){

    // use standard 5 point stencil
    int jp, jm, ip, im;
    int count1=0, count2=0, count3=0, count4=0;

    jp=j-1;
    jm=j+1;

    ip=i-1;
    im=i+1;

    if(ip < 0){
      ip=nx;
      count1+=1;
    }

    if(im > nx){
      im=0;
      count2+=1;
    }

    if(jp < 0){
      jp=ny;
      count3+=1;
    }

    if(jm > ny){
      jm=0;
      count4+=1;
    }

    lap_con[i][j] =(conc[im][j] + conc[ip][j] + conc[i][jm] + conc[i][jp] -4.0*conc[i][j])/(dx*dy);


    // // Test of Periodic Boundary Implementation
    // printf("%d\n", count1 );
    // printf("%d\n", count2 );
    // printf("%d\n", count3 );
    // printf("%d\n", count4 );
    // printf("====\n");
    return lap_con;

  }
