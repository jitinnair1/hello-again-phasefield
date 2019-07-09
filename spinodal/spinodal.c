#include "spinodal.h"

int main(){

  //system Configuration
  int Nx=64,
  Ny=64,
  Nz=0;

  float dx=0.5,
  dy=0.5,
  dz=0.0;

  //time integration parameters:
  int nstep=500,
  istep=50,
  iprint=0;

  float dt=0.01;

  //material Specific Parameters
  float c0 = 0.40,
  noise=0.8,
  mobility = 1.0,
  grad_coef= 0.5;

  float conc[Nx][Ny];
  memset(conc, 0, sizeof(conc[0][0]) * Nx * Ny);

  // *** Check why printing values post memset give non-zero values? ***
  for (size_t i = 0; i <= Nx; i++) {
    for (size_t j = 0; j <= Ny; j++) {
      printf("%f\n", conc[i][j]);
      //conc[i][j]=0.0;
    }
  }
  //get array of random numbers between 0 and 1 for setting initial microstructure
  float random_ZeroToOne_array[Nx][Ny];
  memset(random_ZeroToOne_array, 0, sizeof(random_ZeroToOne_array[0][0]) * Nx * Ny);
  rand_ZeroToOne(Nx, Ny, random_ZeroToOne_array);

  //generate initial microstructure
  prep_microstructure(Nx, Ny, conc, c0, noise, random_ZeroToOne_array);

  printf("Passed\n");
  //write initial microstructure
  write_to_VTK(Nx, Ny, Nz, dx, dy, dz, iprint, conc);

  //evolution loop
  for (iprint=istep; iprint<=nstep; iprint+=istep){

    //get laplacian1
    float lap_con[Nx][Ny];
    //memset(lap_con, 0, sizeof(lap_con[0][0]) * Nx * Ny);
    for (size_t i = 0; i <= Nx; i++) {
      for (size_t j = 0; j <= Ny; j++) {
        lap_con[i][j]=0.0;
      }
    }

    laplacian(Nx, Ny, Nz, dx, dy, dz, conc, lap_con);

    //get Free Energy
    float dfdcon[Nx][Ny];
    //memset(dfdcon, 0, sizeof(dfdcon[0][0]) * Nx * Ny);
    for (size_t i = 0; i <= Nx; i++) {
      for (size_t j = 0; j <= Ny; j++) {
        dfdcon[i][j]=0.0;
      }
    }

    free_energy(Nx, Ny, conc, dfdcon);

    //solve1
    float lap_dummy[Nx][Ny];
    //memset(lap_dummy, 0, sizeof(lap_dummy[0][0]) * Nx * Ny);
    for (size_t i = 0; i <= Nx; i++) {
      for (size_t j = 0; j <= Ny; j++) {
        lap_dummy[i][j]=0.0;
      }
    }

    solve(Nx, Ny, grad_coef, dfdcon, lap_con, lap_dummy);

    //get laplacian2
    float lap2_con[Nx][Ny];
    //memset(lap2_con, 0, sizeof(lap2_con[0][0]) * Nx * Ny);
    for (size_t i = 0; i <= Nx; i++) {
      for (size_t j = 0; j <= Ny; j++) {
        lap2_con[i][j]=0.0;
      }
    }

    laplacian(Nx, Ny, Nz, dx, dy, dz, lap_dummy, lap2_con);

    //solve2
    solve2(Nx, Ny, dt, mobility, lap2_con, conc);

    //write solution to file
    write_to_VTK(Nx, Ny, Nz, dx, dy, dz, iprint, conc);

  }


}
