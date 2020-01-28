#include <stdio.h>
#include <stdlib.h>
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

void write_init_conc(int Nx, int Ny, int NxNy, double conc_print[NxNy]);

void prep_fft(int Nx, int Ny, double dx, double dy, double *kx, double *ky);

double *green_tensor(int Nx, int Ny, double *kx, double *ky, double cm11, double cm12,
                     double cm44, double cp11, double cp12, double cp44, double* tmatx[Nx][Ny][2][2][2][2]);

void elasticity_derivative(int Nx, int Ny, double tmatx, double kx, double ky,
                           double s11, double s22, double s12,double e11, double e22, double e12,
                           double ed11 ,double ed22 ,double ed12 ,double cm11 ,double cm12 ,double cm44 ,
                           double cp11 ,double cp12 ,double cp44 ,double ea, double ei0, fftw_complex* conc) ;
