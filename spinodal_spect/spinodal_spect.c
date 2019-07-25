#include "spinodal_spect.h"

int main(int argc, char const *argv[]) {

  // declarations
  int Nx=128, Ny=128, Nz=0;
  float dx=1.0, dy=1.0, dz=0.0;
  double *kx,*ky;
  fftw_complex *c,*ctilde,*g,*gtilde;
  fftw_plan p1,p2,p3;

  float c0=0.4,
  dt=0.1,
  diffusivity=1.0,
  kappa=1.0,
  A=1.0,
  noise=0.02;

  int nstep=10000,
  iprint=400,
  istep=0;

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
  int ii=0;
  int i, j;

  // get duplicate array and initialise it
  int NxNy=Nx*Ny;
  float conc[NxNy];
  for (size_t i = 0; i < NxNy; i++) {
    conc[i]=0.0;
  }

  //write real part to conc[]
  for (i=0; i < Nx; i++){
    for (j=0; j < Ny; j++){
      ii=i*Nx+j;
      c[ii] = c0 + noise*(0.5-random_ZeroToOne_array[i][j]);
      conc[ii] = creal(c[ii]);
    }
  }

  // write initial concentration to file
  write_to_VTK(Nx, Ny, Nz, dx, dy, dz, istep, NxNy, conc );

  //Boundary condition
  for(i=0; i<Nx; i++){
    if(i < Nx/2)
    kx[i]=2*M_PI*(double)(i)/(double)(Nx*dx);
    else
    kx[i]=2*M_PI*(double)(i-Nx)/(double)(Nx*dx);
  }

  for(j=0; j<Ny; j++){
    if(j<Ny/2)
    ky[j]=2*M_PI*(double)(j)/(double)(Ny*dy);
    else
    ky[j]=2*M_PI*(double)(j-Ny)/(double)(Ny*dy);
  }

  //time loop
  for(istep=1; istep<=nstep; istep++){

    // calculate g
    for(i=0; i<Nx; i++){
      for(j=0; j<Ny; j++){
        ii=i*Nx+j;
        g[ii]=2*A*creal(c[ii])*(1-creal(c[ii]))*(1-2*creal(c[ii]));
      }
    }

    fftw_execute(p3);    //calculating gtilde
    fftw_execute(p1);    //calculating ctilde

    // calculation in fourier space
    for(i=0; i<Nx; i++) {
      for(j=0; j<Ny; j++) {
        ii=i*Nx+j;
        ctilde[ii]=(ctilde[ii]-((kx[i]*kx[i]+ky[j]*ky[j])*
        dt*diffusivity*gtilde[ii]))/
        (1+2*kappa*diffusivity*(kx[i]*kx[i]+ky[j]*ky[j])*
        (kx[i]*kx[i]+ky[j]*ky[j])*dt)
        + _Complex_I*0.0;
      }
    }

    //coming to real space
    fftw_execute(p2);


    //time interval for saving the data
    if(istep%iprint==0)
    {
      //write real part to conc[]
      for (i=0; i < Nx; i++){
        for (j=0; j < Ny; j++){
          ii=i*Nx+j;
          conc[ii] = creal(c[ii])/(double)(Nx*Ny);
        }
      }

      // write initial concentration to file
      write_to_VTK(Nx, Ny, Nz, dx, dy, dz, istep, NxNy, conc );
    }

    //exchanging the values

    for(i=0; i<Nx; i++) {
      for(j=0; j<Ny; j++)
      {
        ii=i*Nx+j;
        c[ii]=creal(c[ii])/(double)(Nx*Ny);
      }
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
  array_deallocate(Ny, random_ZeroToOne_array);

  return 0;
}
