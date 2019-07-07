void write_to_VTK(nx, ny, nz, conc){
  // write data for every istep in VTK file format
  int Nx=nx+1,
  Ny=ny+1,
  Nz=nz+1,
  num_points=0;

  num_points=Nx*Ny*Nz;

  // create output directory if not created
  struct stat st = {0};
  if (stat("./output", &st) == -1) {
    mkdir("./output", 0700);
  }

  FILE *fp;
  char filename[30];
  sprintf(filename, "./output/time_%05d.vtk", k);
  fp = fopen(filename,"w");

  // write header of VTK file
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "time_10.vtk\n");
  fprintf(fp, "ASCII\n");

  // co-ordinates of grid points
  fprintf(fp, "DATASET STRUCTURED_GRID\n");
  fprintf(fp, "DIMENSIONS %6d %6d %6d\n", Nx, Ny, Nz );
  fprintf(fp, "POINTS %6d  float\n", num_points );
  int p, q, r;
  for (p=0; p< Nx; p++){
    for (q=0; q< Ny; q++){
      for (r=0; r< Nz; r++){
        float x, y, z;
        x = p*dx;
        y = q*dy;
        z = r*dz;
        fprintf(fp, "%f %f %f\n", x, y, z);
      }
    }
  }
  // grid point values
  fprintf(fp,"POINT_DATA %6d\n", num_points);
  fprintf(fp,"SCALARS CONC  float  1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for (p=0; p < Nx; p++) {
    fprintf(fp, "%f\n", conc[p]);
  }

  fclose(fp);
}