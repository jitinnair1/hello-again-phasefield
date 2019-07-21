#include "spinodal_spect.h"

int main(int argc, char const *argv[]) {

  // declarations
  int Nx=128, Ny=128;
  double *kx,*ky;
  fftw_complex *c,*ctilde,*g,*gtilde;
  fftw_plan p1,p2,p3;

  float c0=0.4,
  noise=0.02;

  // FFTW allocations
  c=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  ctilde=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));

  g=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  gtilde= (fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));

  kx=(double*)malloc(sizeof(double)*Nx);
  ky=(double*)malloc(sizeof(double)*Ny);

  //creating plans
  p1=fftw_plan_dft_2d(Nx,Ny,c,ctilde,FFTW_FORWARD,FFTW_ESTIMATE);
  p2=fftw_plan_dft_2d(Nx,Ny,ctilde,c,FFTW_BACKWARD,FFTW_ESTIMATE);
  p3=fftw_plan_dft_2d(Nx,Ny,g,gtilde,FFTW_FORWARD,FFTW_ESTIMATE);

  //allocate random array
  float **random_ZeroToOne_array=0;
  random_ZeroToOne_array=array_allocate(Nx, Ny, random_ZeroToOne_array);

  // get array of random numbers between 0 and 1 for setting initial microstructure
  rand_ZeroToOne(Nx, Ny, random_ZeroToOne_array);

  //initial concentrarion
  prep_microstructure(Nx, Ny, c, c0, noise, random_ZeroToOne_array);

  // print conc values to check if random values are generated
  int ii;
  for (size_t i = 0; i < Nx; i++) {
    for (size_t j = 0; j < Ny; j++) {
      ii = i*Nx + j;
      c[ii]=creal(c[ii])/(double)(Nx*Ny);
      printf("%e\n", c[ii] );
    }
  }



  //FFTW Free, Destroy and Cleanup
  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
  fftw_destroy_plan(p3);
  fftw_free(c);
  fftw_free(ctilde);
  fftw_free(g);
  fftw_free(gtilde);
  fftw_cleanup();
  free(kx);
  free(ky);

  return 0;
}
