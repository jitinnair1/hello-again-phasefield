#include "spinodal.h"

int main(){

  //system Configuration
  int Nx=199,
  Ny=199,
  Nz=0;

  float dx=1.0,
  dy=1.0,
  dz=0.0;

  //time integration parameters:
  int nstep=10000,
  iprint=1000,
  istep=0;

  float dt=0.01;

  //material specific parameters
  float c0 = 0.40,
  noise=0.02,
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

  // ***TEST*** : Check size of all arrays and if they are initialized to zero


  //get array of random numbers between 0 and 1 for setting initial microstructure
  rand_ZeroToOne(Nx, Ny, random_ZeroToOne_array);


  // ***TEST*** : Check if Nx * Ny random numbers are being generated between 0 and 1


  //generate initial microstructure
  prep_microstructure(Nx, Ny, conc, c0, noise, random_ZeroToOne_array);

  //write initial microstructure
  write_to_VTK(Nx, Ny, Nz, dx, dy, dz, istep, conc);

  //evolution loop
  for (istep=1; istep<=nstep; istep+=1){

    for (int i=0; i<=Nx; i++){
      for (int j = 0; j <=Ny; j++) {

        //get laplacian1
        lap_con=laplacian(Nx, Ny, Nz, dx, dy, dz, conc, lap_con, i, j);

        // ***TEST*** : Check if PBC is implemented correctly for correct array size

        //get Free Energy
        dfdcon=free_energy(Nx, Ny, conc, dfdcon, i, j);

        // ***TEST*** : Check if free enrgy is computed correctly for correct array size

        //solve1
        lap_dummy=solve(Nx, Ny, grad_coef, dfdcon, lap_con, lap_dummy, i, j);

        // ***TEST*** : Check if constant values make sense and coputaion is done correctly

        //get laplacian2
        lap2_con=laplacian(Nx, Ny, Nz, dx, dy, dz, lap_dummy, lap2_con, i, j);

        //solve2
        conc=solve2(Nx, Ny, dt, mobility, lap2_con, conc, i, j);
        
      }
    }

    if(istep % iprint == 0){
      //write solution to file
      write_to_VTK(Nx, Ny, Nz, dx, dy, dz, istep, conc);
    }

  }

  // deallocate memory
  array_deallocate(Ny, conc);
  array_deallocate(Ny, random_ZeroToOne_array);
  array_deallocate(Ny, lap_dummy);
  array_deallocate(Ny, lap_con);
  array_deallocate(Ny, lap2_con);
  array_deallocate(Ny, dfdcon);

  // ***TEST*** : Check if deallocation is working


  return 0;
}
