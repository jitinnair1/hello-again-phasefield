#include "stdio.h"
#include "math.h"

int main() {

  //Declaring variables
  double alpha, D, dx, dt;
  int j, k, i, l, istep, N, nstep;
  FILE *fp;

  //Initial system configuration
  N=127;
  D=1.0;
  dx=0.5;
  dt=0.1;
  istep=1500;
  nstep=12000;
  alpha=D*dt/(dx*dx);

  //Initial conditions and profile
  float conc[N];
  conc[0]=1.0;
    for (i=1; i<=N; i++){
        conc[i]=0.0;
    }

  //Evolution of profile

  //loop for saving results
  for (k=0; k<=nstep; k+=istep) {

    //loop for evolution
    for (j=0; j<=istep; j++) {

      //loop for explicit condition
      for (i = 0; i < N; i++) {
        conc[i]=(1-2*alpha)*conc[i]+alpha*(conc[i-1]+conc[i+1]);
      }

      //handling boundary condition at right end
      conc[N]=(1-2*alpha)*conc[N]+2*alpha*(conc[N-1]);
    }

      //Write results to file in fpointer after istep iterations
      int Nx=N,
          Ny=1,
          Nz=1,
          num_points=Nx*Ny*Nx;

      char filename[30];
      sprintf(filename, "time_%05d.vtk", k);
      fp = fopen(filename,"w");

      // write header of VTK file
      fprintf(fp, "# vtk DataFile Version 2.0\n");
      fprintf(fp, "time_10.vtk\n");
      fprintf(fp, "ASCII\n");
      fprintf(fp, "DATASET STRUCTURED_GRID\n");

      // co-ordinates of grid points
      fprintf(fp, "DIMENSIONS %6d %6d %6d\n", Nx, Ny, Nz );
      fprintf(fp, "POINTS   %6d  float\n", num_points );
      for (l=0; l<N; l++){
        float x,
              y=0.0,
              z=0.0;
        x = l*dx;
        fprintf(fp, "%f %f %f\n", x, y, z);
      }

      // grid point values
      fprintf(fp,"POINT_DATA %5d\n", num_points);
      fprintf(fp,"SCALARS CONC  float  1\n");
      fprintf(fp,"LOOKUP_TABLE default\n");
      for (l=0; l < N; l++) {
        fprintf(fp, "%f\n", conc[l]);
      }

      fclose(fp);

  }

}
