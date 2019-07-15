#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

// Funcion Declaration
void free_energy(int nx, int ny,
  float conc[nx][ny], float dfdcon[nx][ny]);

void laplacian(
  int nx, int ny, int nz,
  float dx, float dy, float dz,
  float conc[nx][ny], float lap_con[nx][ny]);

void prep_microstructure(int nx, int ny, float conc[nx][ny], float c0,
  float noise, float random_ZeroToOne_array[nx][ny] );

void rand_ZeroToOne(int nx, int ny, float random_ZeroToOne_array[nx][ny]);

void set_array_zero(int nx, int ny, float input_array[nx][ny]);

void solve(int nx, int ny, float grad_coef,
  float dfdcon[nx][ny], float lap_con[nx][ny], float lap2_con[nx][ny]);

void solve2(int nx, int ny, float dt,
  float mobility, float lap2_con[nx][ny], float conc[nx][ny]);

void write_to_VTK( int nx, int ny, int nz,
  float dx, float dy, float dz,
  int iprint, float conc[nx][ny]);
