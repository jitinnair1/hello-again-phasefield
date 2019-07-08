#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>

// Funcion Declaration
void free_energy(int nx, int ny,
  float conc[nx][ny], float dfdcon[nx][ny]);

void laplacian(
  int nx, int ny, int nz,
  float dx, float dy, float dz,
  float conc[nx][ny], float lap_con[nx][ny]);

double rand_ZeroToOne();

void solve(int nx, int ny, float grad_coef,
  float dfdcon[nx][ny], float lap_con[nx][ny], float lap2_con[nx][ny]);

void solve2(int nx, int ny, float dt,
  float mobility, float lap2_con[nx][ny], float conc[nx][ny]);

void write_to_VTK( int nx, int ny, int nz,
  float dx, float dy, float dz,
  int iprint, float conc[nx][ny]);
