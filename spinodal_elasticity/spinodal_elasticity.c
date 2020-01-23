#include "spinodal_elasticity.h"

int main(int argc, char const *argv[]) {

    //start clock
    clock_t tic = clock();

    //declarations
    int Nx=128, Ny=128, Nz=0;
    double dx=1.0, dy=1.0, dz=0.0;
    double *kx,*ky,
            *conc_print, *random_ZeroToOne_array;
    fftw_complex *conc,*conc_tilde,*free_energy,*free_energy_tilde;
    fftw_plan p1,p2,p3;

    double conc0=0.5,
            dt=0.1,
            diffusivity=1.0,
            kappa=1.0,
            A=1.0,
            noise=0.02;

    int nstep=1000,
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

    kx=(double*)malloc(sizeof(double)*Nx);
    ky=(double*)malloc(sizeof(double)*Ny);

    conc_print=(double*)malloc(sizeof(double)*NxNy);
    random_ZeroToOne_array=(double*)malloc(sizeof(double)*NxNy);

    //creating plans
    p1=fftw_plan_dft_2d(Nx, Ny, conc, conc_tilde, FFTW_FORWARD, FFTW_ESTIMATE);
    p2=fftw_plan_dft_2d(Nx, Ny, conc_tilde, conc, FFTW_BACKWARD, FFTW_ESTIMATE);
    p3=fftw_plan_dft_2d(Nx, Ny, free_energy, free_energy_tilde, FFTW_FORWARD, FFTW_ESTIMATE);

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

    //Periodic boundary conditions
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

        //green tensor
        for(i=0; i<Nx; i++){
            for(j=0; j<Ny; j++){
                ii=i*Nx+j;
                free_energy[ii]= 2 * A * creal(conc[ii]) * (1 - creal(conc[ii])) * (1 - 2 * creal(conc[ii]));
            }
        }

        //calculate derivative of elastic_energy
        for(i=0; i<Nx; i++){
            for(j=0; j<Ny; j++){
                ii=i*Nx+j;
                free_energy[ii]= 2 * A * creal(conc[ii]) * (1 - creal(conc[ii])) * (1 - 2 * creal(conc[ii]));
            }
        }

        //calculate derivative of free_energy
        for(i=0; i<Nx; i++){
            for(j=0; j<Ny; j++){
                ii=i*Nx+j;
                free_energy[ii]= 2 * A * creal(conc[ii]) * (1 - creal(conc[ii])) * (1 - 2 * creal(conc[ii]));
            }
        }

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
