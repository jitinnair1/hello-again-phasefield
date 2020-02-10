#include "spinodal_elasticity.h"

void prep_microstructure(int iflag, int Nx, int Ny, int NxNy, fftw_complex conc[NxNy], double conc_print[NxNy], double c0,
                         double noise, double random_ZeroToOne_array[NxNy]) {
  int ii=0;
  int i, j,
  radius=Nx/6,
  cx=Nx/2,
  cy=Ny/2;
  double dist;

  // Initialise with random fluctuation about conc0
  if(iflag==1){
      for (i=0; i < Nx; i++){
          for (j=0; j < Ny; j++){
              ii=i*Nx+j;
              dist = sqrt((i-cx)*(i-cx) + (j-cy)*(j-cy));
              if (dist < radius){
                  conc[ii]=1.0;
              }
              conc_print[ii] = creal(conc[ii]);
          }
      }
  }

  // For running 2D code in 1D mode
  if(iflag==2){
      for (i=0; i < Nx; i++){
          for (j=0; j < Ny; j++){
              ii= i * Nx + j;
              conc[ii] = c0 + noise*(0.5-random_ZeroToOne_array[ii]);
              conc_print[ii] = creal(conc[ii]);
          }
      }
  }
}
