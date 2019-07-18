#include "spinodal.h"

int main(){

  //system Configuration
  int Nx=99,
  Ny=99,
  Nz=0;

  float dx=0.5,
  dy=0.5,
  dz=0.0;

  //time integration parameters:
  int nstep=1000,
  istep=50,
  iprint=0;

  float dt=0.001;

  //material Specific Parameters
  float c0 = 0.40,
  noise=0.8,
  mobility = 1.0,
  grad_coef= 0.5;

  // float **conc;
  // conc=array_allocate(Nx, Ny, conc);
  //
  // float **random_ZeroToOne_array;
  // random_ZeroToOne_array=array_allocate(Nx, Ny, random_ZeroToOne_array);

  int m=Nx+2;
  float **conc = calloc(m, sizeof(float *));
  int r;
  for (r = 0; r <= Ny; r++)
  conc[r] = calloc(Ny, sizeof(float));

  float **random_ZeroToOne_array = calloc(m, sizeof(float *));
  for (r = 0; r <= Ny; r++)
  random_ZeroToOne_array[r] = calloc(Ny, sizeof(float));

  float **lap_dummy = calloc(m, sizeof(float *));
  for (r = 0; r <= Ny; r++)
  lap_dummy[r] = calloc(Ny, sizeof(float));

  float lap_con[Nx][Ny];
  for (size_t i = 0; i < Nx; i++) {
    for (size_t j = 0; j < Ny; j++) {
      lap_con[i][j]=0.0;
    }
  }

  float dfdcon[Nx][Ny];
    for (size_t i = 0; i < Nx; i++) {
      for (size_t j = 0; j < Ny; j++) {
        dfdcon[i][j]=0.0;
      }
    }

    float lap2_con[Nx][Ny];
    for (size_t i = 0; i < Nx; i++) {
      for (size_t j = 0; j < Ny; j++) {
        lap2_con[i][j]=0.0;
      }
    }



  //get array of random numbers between 0 and 1 for setting initial microstructure
  rand_ZeroToOne(Nx, Ny, random_ZeroToOne_array);

  //generate initial microstructure
  prep_microstructure(Nx, Ny, conc, c0, noise, random_ZeroToOne_array);

  //write initial microstructure
  write_to_VTK(Nx, Ny, Nz, dx, dy, dz, iprint, conc);


  //evolution loop
  for (iprint=istep; iprint<=nstep; iprint+=istep){

    //get laplacian1

    laplacian(Nx, Ny, Nz, dx, dy, dz, conc, lap_con);

    //get Free Energy
    free_energy(Nx, Ny, conc, dfdcon);

    //solve1
    solve(Nx, Ny, grad_coef, dfdcon, lap_con, lap_dummy);

    //get laplacian2

    laplacian(Nx, Ny, Nz, dx, dy, dz, lap_dummy, lap2_con);

    //solve2
    solve2(Nx, Ny, dt, mobility, lap2_con, conc);

    //write solution to file
    write_to_VTK(Nx, Ny, Nz, dx, dy, dz, iprint, conc);



  }

  // deallocate memory
  array_deallocate(Ny, conc);
  array_deallocate(Ny, random_ZeroToOne_array);
  array_deallocate(Ny, lap_dummy);


  return 0;
}
