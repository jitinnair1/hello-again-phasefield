#include "spinodal.h"

int main(){

  //System Configuration
  int Nx=100,
  Ny=100,
  Nz=0,
  NxNy=Nx*Ny;

  float dx=0.5,
  dy=0.5,
  dz=0.0;

  // Time integration parameters:
  int nstep=500,
  istep=50;
  iprint=0;

  float dt=0.01;

  // Material Specific Parameters
  float c0 = 0.40,
  mobility = 1.0,
  grad_coef= 0.5,
  conc[Nx][Ny]={0.0},



  // generate initial microstructure
  float noise=0.02;
  int i, j;
  for (i=0; i<=Nx; i++){
    for (j=0; j<=Ny; j++){
      conc[i][j] =c0 + noise*( 0.5-rand_ZeroToOne() );
    } // end for Nx
  } // end for Ny

  //write initial microstructure
  write_to_VTK(Nx, Ny, Nz, dx, dy, dz, iprint, conc);

  //evolution loop
  for (iprint=istep; iprint<=nstep; iprint+=istep){


    //get laplacian1
    lap_con[Nx][Ny]={0.0};
    laplacian(Nx, Ny, Nz, dx, dy, dz, conc, lap_con);

    //get Free Energy
    float dfdcon[Nx][Ny] = {0.0}
    free_energy(Nx, Ny, conc, dfdcon);

    //solve1
    lap_dummy[Nx][Ny]={0.0};
    solve(Nx, Ny, grad_coef, dfdcon, lap_con, lap_dummy);

    //get laplacian2
    lap2_con[Nx][Ny]={0.0};
    laplacian(Nx, Ny, Nz, dx, dy, dz, lap_dummy, lap2_con);

    //solve2
    solve2(Nx, Ny, dt, mobility, lap2_con, conc);

    //write solution to file
    write_to_VTK(nx, ny, nz, dx, dy, dz, iprint, conc);

  }


}
