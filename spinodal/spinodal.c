#include "stdio.h"
#include "math.h"

// Generate random number between 0 and 1
double rand_ZeroToOne()
{
  return (double)rand() / (double)RAND_MAX ;
}

// Generate discrete laplacian
double laplacian(Nx, Ny, Nz, dx, dy, dz, &conc){

  // use standard 5 point stencil

}

// main function

int main() {
  //Declararions
  int j, k, i, istep, N;
  FILE *fp;

  //System Configuration
  int Nx=100,
  Ny=100,
  NxNy=Nx*Ny;

  float dx=0.5,
  dy=0.5;

  // Time integration parameters:
  int nstep=50,
  nprint=50;

  float dtime=1.0E-4;

  // Material Specific Parameters
  float c0 = 0.40,
  mobility = 1.0,
  grad_coef= 0.5,
  conc[Nx][Ny]={0.0};

  // generate initial microstructure
  float noise=0.02;

  for (i=0, i<Nx; i++)  {
    for (i=0, i<Ny; i++) {
      conc[i][j] =c0 + noise*(0.5-rand_ZeroToOne());
    } // end for Nx
  } // end for Ny


  // get laplacian

  // get Free Energy

  // solve

  // write solution to file

  //
  for(i = 0; i < nstep; i+=nprint) {
    char filename[30];
    sprintf(filename, "file%06d.txt", i);
    fp = fopen(filename,"w");

    // write array in loop after creating VTK template

    fclose(fp);
  }


}
