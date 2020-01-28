//
// Created by Jitin Nair on 24/01/20.
//

#include "spinodal_elasticity.h"


void elasticity_derivative(int Nx, int Ny, double tmatx, double kx, double ky,
                           double s11, double s22, double s12,double e11, double e22, double e12,
                           double ed11 ,double ed22 ,double ed12 ,double cm11 ,double cm12 ,double cm44 ,
                           double cp11 ,double cp12 ,double cp44 ,double ea, double ei0, fftw_complex* conc,
                           fftw_plan p4, fftw_plan p10 ) {

    int niter=10, NxNy=Nx*Ny;
    double tolerance=0.001;

    for (int ii = 0; ii < NxNy; ++ii) {
        ei11[ii] = ei0*conc[ii];
        ei22[ii] = ei0*conc[ii];
        ei33[ii] = ei0*conc[ii];
        ei12[ii] = 0.0*conc[ii];

        c11[ii] = conc[ii]*cp11 +(1.0-conc[ii])*cm11;
        c12[ii] = conc[ii]*cp12 +(1.0-conc[ii])*cm12;
        c44[ii] = conc[ii]*cp44 +(1.0-conc[ii])*cm44;
    }

    for (int k = 0; k < niter; ++k) {

        //take stress and strains to fourier space
        fftw_execute_dft(p4, e11, e11k);
        fftw_execute_dft(p4, e12, e12k);
        fftw_execute_dft(p4, e22, e22k);
        fftw_execute_dft(p4, s11, s11k);
        fftw_execute_dft(p4, s12, s12k);
        fftw_execute_dft(p4, s22, s22k);

        //assemble ematx and smatx
        for (int i = 0; i < Nx; ++i) {

        }

        //get updated values of strains in fourier space

        //take strains to real space

        //calculate stress in real space

        //check for convergence


    }

    //return value of derivative

}

