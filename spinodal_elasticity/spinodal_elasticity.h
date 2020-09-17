#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gperftools/profiler.h>

void prep_microstructure(int iflag, int Nx, int Ny, int NxNy, fftw_complex conc[NxNy], double conc_print[NxNy], double c0,
                         double noise, double random_ZeroToOne_array[NxNy] ) ;

double *rand_ZeroToOne(int Nx, int Ny, int seed, double *random_ZeroToOne_array) ;

void write_to_VTK( int nx, int ny, int nz,
                   double dx, double dy, double dz,
                   int iprint, int NxNy, double conc[NxNy] );

void write_init_conc(int Nx, int Ny, int NxNy, double conc_print[NxNy], char label[30]);

void prep_fft(int Nx, int Ny, double dx, double dy, double *kx, double *ky,
              double *k2, double *k4);

double *green_tensor(int Nx, int Ny, double kx[], double ky[], double cm11, double cm12,
                     double cm44, double cp11, double cp12, double cp44, double tmatx[Nx][Ny][2][2][2][2]);

void elasticity_derivative(int Nx, int Ny, int num_points,
                           double tmatx[][Ny][2][2][2][2], double smatx[][Ny][2][2], double ematx[][Ny][2][2],
                           fftw_complex s11[num_points], fftw_complex s22[num_points], fftw_complex s12[num_points],
                           fftw_complex e11[num_points], fftw_complex e22[num_points], fftw_complex e12[num_points],
                           fftw_complex s11k[num_points], fftw_complex s22k[num_points], fftw_complex s12k[num_points],
                           fftw_complex e11k[num_points], fftw_complex e22k[num_points], fftw_complex e12k[num_points],
                           double ed11[num_points], double ed22[num_points], double ed12[num_points],
                           double et11[num_points], double et22[num_points], double et12[num_points],
                           double ei11[num_points], double ei22[num_points], double ei33[num_points], double ei12[num_points],
                           double cm11, double cm12, double cm44,
                           double c11[num_points], double c12[num_points], double c44[num_points],
                           double cp11, double cp12, double cp44, double ea[], double ei0,
                           fftw_complex conc[num_points], fftw_complex delsdc[num_points], fftw_plan p4, fftw_plan p5,
                           int istep);