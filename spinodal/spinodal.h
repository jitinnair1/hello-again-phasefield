#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

// Funcion Declaration
void free_energy(int nx, int ny, float conc[][], float dfdcon[][]);

void laplacian(
  int nx, int ny, int nz,
  float dx, float dy, float dz,
  float conc[][], int lap_con[][]);

double rand_ZeroToOne();

void solve(int nx, int ny, float grad_coef,
  float dfdcon[][], float lap_con[][], lap2_con[][]);

void solve2(int nx, int ny, float dt,
  float mobility, float lap2_con, float conc);

void write_to_VTK( int nx, int ny, int nz,
  float dx, float dy, float dz,
  int iprint, float conc[][]);
