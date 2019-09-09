#include "ternary_coarsening.h"

int main(int argc, char const *argv[]) {

  //start clock
  clock_t tic = clock();

  //declarations for FFTW
  double *kx,*ky, *k2, *k4,
  *mobility_BB, *mobility_CC, *mobility_BC,
  *concA_print, *concB_print, *concC_print;

  fftw_complex *conc_A, *conc_B, *conc_C,
  *conc_B_tilde, *conc_C_tilde,
  *conc_B_tilde_new, *conc_C_tilde_new,
  *dfree_dB, *dfree_dC,
  *dfree_dB_tilde, *dfree_dC_tilde;

  fftw_plan p1, p2, p3, p4, p5, p6;

  //system specific
  int Nx=512, Ny=512, Nz=0;
  double dx=1.0, dy=1.0, dz=0.0;
  double dt=0.001;
  double conc_C0=0.333,
  kappa_A=2.5,
  kappa_B=2.5,
  kappa_C=2.5,
  kappa_BB=kappa_A+kappa_B,
  kappa_CC=kappa_A+kappa_C,
  kappa_BC=kappa_A,
  kai_AB=5.0,
  kai_AC=5.0,
  kai_BC=5.0,
  mobility_A=1.0,
  mobility_B=1.0,
  mobility_C=1.0,
  noise=0.02;

  //for output printing and evolution-loop control
  int nstep=1,
  iprint=1,
  istep=0;

  //others
  int ii=0;
  int i, j;
  int NxNy=Nx*Ny;
  double detL,
  L1,
  L2,
  L3,
  L4,
  R1,
  R2;

  //FFTW allocations
  conc_A=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  conc_B=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  conc_C=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  conc_B_tilde=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  conc_C_tilde=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));

  conc_B_tilde_new=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  conc_C_tilde_new=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));

  dfree_dB=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  dfree_dB_tilde=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));

  dfree_dC=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  dfree_dC_tilde=(fftw_complex*)fftw_malloc(Nx*Ny*sizeof(fftw_complex));

  mobility_BB=(double*)malloc(Nx*Ny*sizeof(double));
  mobility_CC=(double*)malloc(Nx*Ny*sizeof(double));
  mobility_BC=(double*)malloc(Nx*Ny*sizeof(double));


  kx=(double*)malloc(sizeof(double)* Nx);
  ky=(double*)malloc(sizeof(double)* Ny);
  k2=(double*)malloc(sizeof(double)* NxNy);
  k4=(double*)malloc(sizeof(double)* NxNy);

  //creating plans
  p1=fftw_plan_dft_2d(Nx,Ny,conc_B,conc_B_tilde,FFTW_FORWARD,FFTW_ESTIMATE);
  p2=fftw_plan_dft_2d(Nx,Ny,conc_C,conc_C_tilde,FFTW_FORWARD,FFTW_ESTIMATE);
  p3=fftw_plan_dft_2d(Nx,Ny,conc_B_tilde,conc_B,FFTW_BACKWARD,FFTW_ESTIMATE);
  p4=fftw_plan_dft_2d(Nx,Ny,conc_C_tilde,conc_C,FFTW_BACKWARD,FFTW_ESTIMATE);
  p5=fftw_plan_dft_2d(Nx,Ny,dfree_dB,dfree_dB_tilde,FFTW_FORWARD,FFTW_ESTIMATE);
  p6=fftw_plan_dft_2d(Nx,Ny,dfree_dC,dfree_dC_tilde,FFTW_FORWARD,FFTW_ESTIMATE);

  //allocate random array
  float **random_ZeroToOne_array1=0,
  **random_ZeroToOne_array2=0;
  random_ZeroToOne_array1=array_allocate(Nx, Ny, random_ZeroToOne_array1);
  random_ZeroToOne_array2=array_allocate(Nx, Ny, random_ZeroToOne_array2);

  // get array of random numbers between 0 and 1 for setting initial microstructure
  random_ZeroToOne_array1=rand_ZeroToOne(Nx, Ny, random_ZeroToOne_array1);
  random_ZeroToOne_array2=rand_ZeroToOne(Nx, Ny, random_ZeroToOne_array2);

  // get duplicate array and initialise it
  concA_print = (double *) malloc(sizeof(double) * NxNy);
  concB_print = (double *) malloc(sizeof(double) * NxNy);
  concC_print = (double *) malloc(sizeof(double) * NxNy);

  for (i = 0; i < NxNy; i++) {
    concA_print[i]=0.0;
    concB_print[i]=0.0;
    concC_print[i]=0.0;
  }

  //write real part to conc[]
  for (i=0; i < Nx; i++){
    for (j=0; j < Ny; j++){
      ii=i*Nx+j;
      conc_C[ii] = conc_C0 + noise*(0.5-random_ZeroToOne_array1[i][j]);
      conc_B[ii] = conc_C[ii];
      conc_A[ii] = 1-conc_B[ii]-conc_C[ii];
      concA_print[ii] = creal(conc_A[ii]);
      concB_print[ii] = creal(conc_B[ii]);
      concC_print[ii] = creal(conc_C[ii]);
    }
  }

  // write initial concentration to file
  write_to_VTK(Nx, Ny, Nz, dx, dy, dz, istep, NxNy, concA_print, concB_print, concC_print );

  //print completion status
  printf("Timestep %d completed\n", istep );


  //time loop
  for(istep=1; istep<=nstep; istep++){

    // calculate g
    for(i=0; i<Nx; i++){
      for(j=0; j<Ny; j++){
        ii=i*Nx+j;
        dfree_dB[ii]=log(creal(conc_B[ii]))-log(creal(conc_A[ii]))
        + (kai_BC-kai_AC)*creal(conc_C[ii]) + kai_AB*(creal(conc_A[ii])-creal(conc_B[ii]));
        dfree_dC[ii]=log(creal(conc_C[ii]))-log(creal(conc_A[ii]))
        + (kai_BC-kai_AB)*creal(conc_B[ii]) + kai_AC*(creal(conc_A[ii])-creal(conc_C[ii]));
      }
    }

    fftw_execute(p1);
    fftw_execute(p2);
    fftw_execute(p5);
    fftw_execute(p6);

    // calculation in fourier space
    for(i=0; i<Nx; i++) {
      for(j=0; j<Ny; j++) {

        //Boundary condition
        if(i < Nx/2)
            kx[i]=2*M_PI*(double)(i)/(double)(Nx*dx);
        else
            kx[i]=2*M_PI*(double)(i-Nx)/(double)(Nx*dx);

        if(j<Ny/2)
            ky[j]=2*M_PI*(double)(j)/(double)(Ny*dy);
        else
            ky[j]=2*M_PI*(double)(j-Ny)/(double)(Ny*dy);

        ii=i*Nx+j;
        //get k-vectors
        k2[ii]=(kx[i]*kx[i]+ky[j]*ky[j]);
        k4[ii]=k2[ii]*k2[ii];

        //update mobility values
        mobility_BB[ii]=(1-creal(conc_B[ii]))*(1-creal(conc_B[ii]))*mobility_B
                          + creal(conc_B[ii])*creal(conc_B[ii])*(mobility_A+mobility_C);
        mobility_CC[ii]=(1-creal(conc_C[ii]))*(1-creal(conc_C[ii]))*mobility_C
                          + creal(conc_C[ii])*creal(conc_C[ii])*(mobility_A+mobility_B);
        mobility_BC[ii]=creal(conc_B[ii])*creal(conc_C[ii])*mobility_A
                          - (1-creal(conc_C[ii]))*creal(conc_B[ii])*mobility_C
                          + (1-creal(conc_B[ii]))*creal(conc_C[ii])*mobility_B;

        //determinant
        L1=1 + 2*mobility_BB[ii]*kappa_BB*k4[ii]*dt
        - 2*mobility_BC[ii]*kappa_BC*k4[ii]*dt;

        L2=2*(mobility_BB[ii]*kappa_BC - mobility_BC[ii]*kappa_CC)*k4[ii]*dt;

        L3=2*(mobility_CC[ii]*kappa_BC - mobility_BC[ii]*kappa_BB)*k4[ii]*dt;

        L4=1 + 2*mobility_CC[ii]*kappa_CC*k4[ii]*dt
        - 2*mobility_BC[ii]*kappa_BC*k4[ii]*dt;

        detL = L1*L4 - L2*L3;

        //inverted matrix parameters
        R1=conc_B_tilde[ii] - mobility_BB[ii]*k2[ii]*dfree_dB_tilde[ii]*dt
        + mobility_BC[ii]*k2[ii]*dfree_dC[ii]*dt;
        R2=conc_C_tilde[ii] - mobility_CC[ii]*k2[ii]*dfree_dC_tilde[ii]*dt
        + mobility_BC[ii]*k2[ii]*dfree_dB[ii]*dt;

        //update conc_B
        conc_B_tilde_new[ii]=(L4*R1 - R2*L2)/detL;
        conc_B_tilde[ii]=conc_B_tilde_new[ii];

        //update CONC_C
        conc_C_tilde_new[ii]=(-L3*R1 + R2*L1)/detL;
        conc_C_tilde[ii]=conc_C_tilde_new[ii];
      }
    }

    //coming to real space
    fftw_execute(p3);
    fftw_execute(p4);


    //time interval for saving the data
    if(istep%iprint==0)
    {
      //write real part to conc[]
      for (i=0; i < Nx; i++){
        for (j=0; j < Ny; j++){
          ii=i*Nx+j;
          concB_print[ii] = creal(conc_B[ii])/(double) NxNy;
          concC_print[ii] = creal(conc_C[ii])/(double) NxNy;
          concA_print[ii] = 1-concB_print[ii]-concC_print[ii];
        }
      }

      // write initial concentration to file
      write_to_VTK(Nx, Ny, Nz, dx, dy, dz, istep, NxNy, concA_print, concB_print, concC_print );

      //print completion status
      printf("Timestep %d completed\n", istep );
    }

    //set complex part of conc_arrays to zero
    for(i=0; i<Nx; i++) {
      for(j=0; j<Ny; j++)
      {
        ii=i*Nx+j;
        conc_C[ii] = creal(conc_C[ii])/(double) NxNy;
        conc_B[ii] = creal(conc_B[ii])/(double) NxNy;
        conc_A[ii] = 1-conc_C[ii]-conc_B[ii];
      }
    }
  }

  //FFTW Free, Destroy and Cleanup
  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
  fftw_destroy_plan(p3);
  fftw_destroy_plan(p4);
  fftw_destroy_plan(p5);
  fftw_destroy_plan(p6);

  free(mobility_BB);
  free(mobility_CC);
  free(mobility_BC);
  free(k2);free(k4);
  free(kx);free(ky);
  array_deallocate(Ny, random_ZeroToOne_array1);
  array_deallocate(Ny, random_ZeroToOne_array2);
  fftw_free(conc_B_tilde);
  fftw_free(conc_B_tilde_new);
  fftw_free(conc_C_tilde);
  fftw_free(conc_C_tilde_new);
  fftw_free(dfree_dB);
  fftw_free(dfree_dC);
  fftw_free(dfree_dB_tilde);
  fftw_free(dfree_dC_tilde);
  fftw_free(conc_A);
  fftw_free(conc_B);
  fftw_free(conc_C);
  fftw_cleanup();

  // end clock
  clock_t toc = clock();
  printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

  return 0;
}
