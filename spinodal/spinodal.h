#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

// Funcion Declaration

//float* array_allocate(int nx, int ny, float **array);

void array_deallocate(int ny, float **array);

void free_energy(int nx, int ny,
  float **conc, float dfdcon[nx][ny]);

void laplacian(
  int nx, int ny, int nz,
  float dx, float dy, float dz,
  float **conc, float lap_con[nx][ny]);

void prep_microstructure(int nx, int ny, float **conc, float c0,
  float noise, float **random_ZeroToOne_array );

void rand_ZeroToOne(int nx, int ny, float **random_ZeroToOne_array);

void set_array_zero(int nx, int ny, float input_array[nx][ny]);

void solve(int nx, int ny, float grad_coef,
  float dfdcon[nx][ny], float lap_con[nx][ny], float **lap_dummy);

void solve2(int nx, int ny, float dt,
  float mobility, float lap2_con[nx][ny], float **conc);

void write_to_VTK( int nx, int ny, int nz,
  float dx, float dy, float dz,
  int iprint, float **conc);
