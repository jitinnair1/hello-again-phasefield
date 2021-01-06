#include "spinodal_spect.h"

void prep_microstructure(int iflag, int Nx, int Ny, int NxNy, fftw_complex conc[NxNy], double conc_print[NxNy], double c0,
                         double noise, double random_ZeroToOne_array[NxNy] ) {
  int ii=0;
  int i, j;

  // Initialise with random fluctuation about conc0
  if(iflag==1){
      for (i=0; i < Nx; i++){
          for (j=0; j < Ny; j++){
              ii= i * Nx + j;
              conc[ii] = c0 + noise*(0.5-random_ZeroToOne_array[ii]);
              conc_print[ii] = creal(conc[ii]);
          }
      }
  }

  // For running 2D code in 1D mode
  if(iflag==2){
      for (i = 0; i < Nx; ++i) {
          ii=i*Nx;
          conc[ii] = c0 + noise * (0.5 - random_ZeroToOne_array[ii]);
          conc_print[ii] = creal(conc[ii]);
      }

      for (i = 0; i < Nx; i++) {
          for (j = 0; j < Ny; ++j) {
              ii=i*Nx+j;
              conc[ii]=conc[i*Nx];
              conc_print[ii] = creal(conc[ii]);
          }
      }
  }
}
