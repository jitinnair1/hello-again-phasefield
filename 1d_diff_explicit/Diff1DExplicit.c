#include "stdio.h"
#include "math.h"

void write_to_VTK(nx, ny, nz, k, conc){

  //Write results to file in fpointer after istep iterations
  int Nx=nx+1,
      Ny=ny+1,
      Nz=nz+1,
      num_points=0;

  num_points=Nx*Ny*Nz;

  char filename[30];
  sprintf(filename, "time_%05d.vtk", k);
  fp = fopen(filename,"w");

  // write header of VTK file
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "time_10.vtk\n");
  fprintf(fp, "ASCII\n");

  // co-ordinates of grid points
  fprintf(fp, "DATASET STRUCTURED_GRID\n");
  fprintf(fp, "DIMENSIONS %6d %6d %6d\n", Nx, Ny, Nz );
  fprintf(fp, "POINTS %6d  float\n", num_points );
  for (l=0; l<=N; l++){
    float x, y, z;
    x = l*dx;
    y = l*dy;
    z = l*dz;
    fprintf(fp, "%f %f %f\n", x, y, z);
  }

  // grid point values
  fprintf(fp,"POINT_DATA %6d\n", num_points);
  fprintf(fp,"SCALARS CONC  float  1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for (l=0; l <= N; l++) {
    fprintf(fp, "%f\n", conc[l]);
  }

  fclose(fp);
}


int main() {

  //Declaring variables
  double alpha, D, dx, dy, dz, dt;
  int j, k, i, l, istep, N, nstep;
  FILE *fp;

  //Initial system configuration
  N=511;
  D=1.0;
  dx=0.8;
  dy=0,0;
  dz=0.0;
  dt=0.01;
  istep=10;
  nstep=500;
  alpha=D*dt/(dx*dx);

  //Initial conditions and profile
  float conc[N];
  conc[0]=1.0;
    for (i=1; i<N; i++){
        conc[i]=0.0;
    }

  write_to_VTK(N, 0, 0, 0, conc);
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

      // Call VTK write function
      write_to_VTK(N, 0, 0, k, conc);

  }

}
