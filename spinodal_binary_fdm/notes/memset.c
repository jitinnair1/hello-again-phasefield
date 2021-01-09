#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(){
  int Nx=100, Ny=100;
  float array[Nx][Ny];
  memset(array, 0, sizeof(array[0][0]) * Nx * Ny);

  // *** Check why printing values post memset give non-zero values? ***
 for (int i = 0; i <= Nx; i++) {
   for (int j = 0; j <= Ny; j++) {
     printf("conc[%d][%d] = %f\n", i, j, array[i][j]);
     //conc[i][j]=0.0;
   }
 }
 return 0;
}
