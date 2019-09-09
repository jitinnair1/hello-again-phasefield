#include "Diff1DExplicit.h"

int main() {

  //Declaring variables
  double alpha, D, dx, dy, dz, dt;
  int j, i, nx, ny, nz, istep, iprint, N, nstep;

  //Initial system configuration

  N=255;                //Loops go from 0 to N both inclusive
  D=1.0;                //Diffusivity
  dx=0.5;
  dy=0.0;
  dz=0.0;
  dt=0.001;
  istep=100;            //timesteps after which values are written to file
  iprint=0;             //variable to print initial state
  nstep=1000;           //total number of timesteps
  alpha=D*dt/(dx*dx);   //parameter used for computation

  //file write parameters
  nx=N;
  ny=0;
  nz=0;

  //Initial conditions and profile
 double conc[N];
  conc[0]=1.0;
  for (i=1; i<N; i++){
    conc[i]=0.0;
  }

  //write initial concentration
  write_to_VTK( nx, ny, nz, dx, dy, dz, iprint, conc );

  //Evolution of profile

  //loop for saving results
  for (iprint=istep; iprint<=nstep; iprint+=istep) {

    //loop for evolution
    for (j=0; j<=istep; j++) {

      //loop for explicit condition
      for (i = 1; i < N; i++) {
        conc[i]=(1-2*alpha)*conc[i]+alpha*(conc[i-1]+conc[i+1]);
      }

      //handling boundary condition at right end
      conc[N]=(1-2*alpha)*conc[N]+2*alpha*(conc[N-1]);
    }

    //Write results to file in fpointer after istep iterations
    write_to_VTK( nx, ny, nz, dx, dy, dz, iprint, conc );
    printf("Completed Timestep %d\n", iprint);

  }
}
