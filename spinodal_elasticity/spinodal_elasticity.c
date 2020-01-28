#include "spinodal_elasticity.h"

int main(int argc, char const *argv[]) {

    //start clock
    clock_t tic = clock();

    //declarations
    int Nx=128, Ny=128, Nz=0;
    double dx=1.0, dy=1.0, dz=0.0;
    double *kx,*ky, *conc_print, *random_ZeroToOne_array;
    double *tmatx, *ed11, *ed22, *ed12;
    fftw_complex *e11, *e22, *e12, *s11, *s22, *s12,
            *e11k, *e22k, *e12k, *s11k, *s22k, *s12k;
    fftw_complex *conc,*conc_tilde,*free_energy,*free_energy_tilde;
    fftw_plan p1,p2,p3,p4,p5,p6,p7,
    p8,p9,p10,p11,p12,p13,p14,p15;

    double conc0=0.5,
            dt=0.05,
            diffusivity=1.0,
            kappa=1.0,
            A=1.0,
            noise=0.02;

    // elastic constants:
    double cm11, cm12, cm44,
    cp11, cp12, cp44;

    cm11 = 1400.0;
    cm12 = 600.0;
    cm44 = 400.0;
    cp11 = 2.0*cm11;
    cp12 = 2.0*cm12;
    cp44 = 2.0*cm44;

    //eigen strains
    double ei0 = 0.01;

    //applied strains
    double ea[3];
    ea[0]=0.0;
    ea[1]=0.01;
    ea[2]=0.0;

    int nstep=5000,
            iprint=100,
            istep=0;

    int i, j, ii;
    int NxNy=Nx*Ny;
    int iflag=1;

    //FFTW allocations
    conc=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));
    conc_tilde=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));

    free_energy=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));
    free_energy_tilde= (fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));

    //initial stess and strain components
    s11=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));
    s22=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));
    s12=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));
    s11k=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));
    s22k=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));
    s12k=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));

    e11=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));
    e22=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));
    e12=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));
    e11k=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));
    e22k=(fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));
    e12k=((fftw_complex*)fftw_malloc(Nx * Ny * sizeof(fftw_complex));


    //FFT related all
    kx=(double*)malloc(sizeof(double)*Nx);
    ky=(double*)malloc(sizeof(double)*Ny);

    //For RNG and write_to_VTK
    conc_print=(double*)malloc(sizeof(double)*NxNy);
    random_ZeroToOne_array=(double*)malloc(sizeof(double)*NxNy);

    //eigen strains
    ed11=(double*)malloc(sizeof(double)*NxNy);
    ed22=(double*)malloc(sizeof(double)*NxNy);
    ed12=(double*)malloc(sizeof(double)*NxNy);

    //green tensor array
    tmatx = (double*) malloc(sizeof(double)*Nx*Ny*2*2*2*2*2*2);

    //creating plans
    p1=fftw_plan_dft_2d(Nx, Ny, conc, conc_tilde, FFTW_FORWARD, FFTW_ESTIMATE);
    p2=fftw_plan_dft_2d(Nx, Ny, conc_tilde, conc, FFTW_BACKWARD, FFTW_ESTIMATE);
    p3=fftw_plan_dft_2d(Nx, Ny, free_energy, free_energy_tilde, FFTW_FORWARD, FFTW_ESTIMATE);

    //strains forward
    p4=fftw_plan_dft_2d(Nx, Ny, e11, e11k, FFTW_FORWARD, FFTW_ESTIMATE);
    p5=fftw_plan_dft_2d(Nx, Ny, e22, e22k, FFTW_FORWARD, FFTW_ESTIMATE);
    p6=fftw_plan_dft_2d(Nx, Ny, e12, e12k, FFTW_FORWARD, FFTW_ESTIMATE);
    p7=fftw_plan_dft_2d(Nx, Ny, s11, s11k, FFTW_FORWARD, FFTW_ESTIMATE);
    p8=fftw_plan_dft_2d(Nx, Ny, s12, s12k, FFTW_FORWARD, FFTW_ESTIMATE);
    p9=fftw_plan_dft_2d(Nx, Ny, s22, s22k, FFTW_FORWARD, FFTW_ESTIMATE);

    //strains backward
    p10=fftw_plan_dft_2d(Nx, Ny, e11k, e11, FFTW_BACKWARD, FFTW_ESTIMATE);
    p11=fftw_plan_dft_2d(Nx, Ny, e22k, e22, FFTW_BACKWARD, FFTW_ESTIMATE);
    p12=fftw_plan_dft_2d(Nx, Ny, e12k, e12, FFTW_BACKWARD, FFTW_ESTIMATE);
    p13=fftw_plan_dft_2d(Nx, Ny, s11k, s11, FFTW_BACKWARD, FFTW_ESTIMATE);
    p14=fftw_plan_dft_2d(Nx, Ny, s12k, s12, FFTW_BACKWARD, FFTW_ESTIMATE);
    p15=fftw_plan_dft_2d(Nx, Ny, s22k, s22, FFTW_BACKWARD, FFTW_ESTIMATE);

    // get array of random numbers between 0 and 1 for setting initial microstructure
    rand_ZeroToOne(Nx, Ny, 0, random_ZeroToOne_array);

    //prepare microstructure
    prep_microstructure(iflag, Nx, Ny, NxNy, conc, conc_print, conc0, noise, random_ZeroToOne_array);

    // write initial concentration to file
    write_to_VTK(Nx, Ny, Nz, dx, dy, dz, istep, NxNy, conc_print);

    //write initial conc in TXT
    write_init_conc(Nx, Ny, NxNy, conc_print);

    //print completion status
    printf("Timestep %d completed\n", istep );

    //prepare kx and ky
    prep_fft(Nx, Ny, dx, dy, kx, ky);

    //green tensor
    tmatx = green_tensor(Nx, Ny, kx, ky, cm11, cm12, cm44,
                         cp11, cp12, cp44, tmatx);

    //time loop
    for(istep=1; istep<=nstep; istep++){


        //calculate derivative of free_energy
        for(i=0; i<Nx; i++){
            for(j=0; j<Ny; j++){
                ii=i*Nx+j;
                free_energy[ii]= 2 * A * creal(conc[ii]) * (1 - creal(conc[ii])) * (1 - 2 * creal(conc[ii]));
            }
        }

        //calculate derivative of elastic_energy
        elasticity_derivative(Nx,Ny,tmatx,kx,ky,
        s11,s22,s12,e11,e22,e12,
        ed11,ed22,ed12,cm11,cm12,cm44,
        cp11,cp12,cp44,ea,ei0,conc);

        fftw_execute(p3);    //calculating free_energy_tilde
        fftw_execute(p1);    //calculating conc_tilde

        // calculation in fourier space
        for(i=0; i<Nx; i++) {
            for(j=0; j<Ny; j++) {
                ii=i*Nx+j;
                conc_tilde[ii]= (conc_tilde[ii] - ((kx[i] * kx[i] + ky[j] * ky[j]) * dt * diffusivity * free_energy_tilde[ii])) /
                                (1+2*kappa*diffusivity*(kx[i]*kx[i]+ky[j]*ky[j])*(kx[i]*kx[i]+ky[j]*ky[j])*dt)
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
                    conc_print[ii] = creal(conc[ii]) / (double)(Nx * Ny);
                }
            }

            // write initial concentration to file
            write_to_VTK(Nx, Ny, Nz, dx, dy, dz, istep, NxNy, conc_print );

            //print completion status
            printf("Timestep %d completed\n", istep );
        }

        //exchanging the values
        for(i=0; i<Nx; i++) {
            for(j=0; j<Ny; j++)
            {
                ii=i*Nx+j;
                conc[ii]= creal(conc[ii]) / (double)(Nx * Ny);
            }
        }
    }

    //FFTW Free, Destroy and Cleanup
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(p3);
    fftw_free(conc);
    fftw_free(conc_tilde);
    fftw_free(free_energy);
    fftw_free(free_energy_tilde);
    free(kx);
    free(ky);
    free(conc_print);
    free(random_ZeroToOne_array);
    fftw_cleanup();

    // end clock
    clock_t toc = clock();
    printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);


    return 0;
}

