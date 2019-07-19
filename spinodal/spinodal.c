#include "spinodal.h"

int main(){

  //system Configuration
  int Nx=199,
  Ny=199,
  Nz=0;

  float dx=0.5,
  dy=0.5,
  dz=0.0;

  //time integration parameters:
  int nstep=25000,
  istep=1000,
  iprint=0;

  float dt=0.001;

  //material Specific Parameters
  float c0 = 0.50,
  noise=1,
  mobility = 1.0,
  grad_coef= 0.5;


  //memory allocation
  float **conc=0,
  **random_ZeroToOne_array=0,
  **lap_con=0,
  **dfdcon=0,
  **lap2_con=0,
  **lap_dummy=0;

  conc=array_allocate(Nx, Ny, conc);
  random_ZeroToOne_array=array_allocate(Nx, Ny, random_ZeroToOne_array);
  lap_dummy=array_allocate(Nx, Ny, lap_dummy);
  lap_con=array_allocate(Nx, Ny, lap_con);
  dfdcon=array_allocate(Nx, Ny, dfdcon);
  lap2_con=array_allocate(Nx, Ny, lap2_con);

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
  array_deallocate(Ny, lap_con);
  array_deallocate(Ny, lap2_con);
  array_deallocate(Ny, dfdcon);
  return 0;
}
