#include "spinodal_elasticity.h"

void write_var_to_VTK( int nx, int ny, int nz,
  double dx, double dy, double dz,
  int iprint, int NxNy, double conc[NxNy], char label[30]) {

  // write data for every istep in VTK file format
  int Nx=nx,
  Ny=ny,
  Nz=nz+1,
  num_points=0;

  num_points=Nx*Ny*Nz;

  // create output directory if not created
  struct stat st = {0};
  if (stat("./output", &st) == -1) {
        mkdir("./output", 0700);
  }

  char path[30] = "./output/";
  char time[30];
  sprintf(time, "_%05d", iprint);
  char filename[strlen(label) + strlen(path) + strlen(time) + 5];
  strcpy(filename, path);
  strcat(filename, label);
  strcat(filename, time);
  strcat(filename, ".vtk");

  FILE *fp;
  fp = fopen(filename,"w");

  // write header of VTK file
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "time_10.vtk\n");
  fprintf(fp, "ASCII\n");

  // co-ordinates of grid points
  fprintf(fp, "DATASET STRUCTURED_GRID\n");
  fprintf(fp, "DIMENSIONS %6d %6d %6d\n", Nx, Ny, Nz );
  fprintf(fp, "POINTS %7d  float\n", num_points );
  int p, q;
  double x, y, z;

  for (p=0; p < Nx ; p++){
    for (q=0; q < Ny; q++){
        x = p*dx;
        y = q*dy;
        z = 0.0;
        fprintf(fp, "%f %f %f\n", x, y, z);
      }
    }


  // grid point values
  fprintf(fp,"POINT_DATA %6d\n", num_points);
  fprintf(fp,"SCALARS %s float  1\n", label);
  fprintf(fp,"LOOKUP_TABLE default\n");
  int i, j, ii;
  for (i = 0; i < Nx; i++) {
    for (j = 0; j < Ny; j++) {
      ii=i*Nx + j;
      fprintf(fp, "%lf\n", conc[ii] );
    }
  }
  fclose(fp);
}
