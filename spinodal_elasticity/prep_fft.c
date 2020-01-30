//
// Created by Jitin Nair on 24/01/20.
//
#include "spinodal_elasticity.h"

void prep_fft(int Nx, int Ny, double dx, double dy, double *kx, double *ky,
        double *k2, double *k4){
    //Periodic boundary conditions
    for(int i=0; i<Nx; i++){
        if(i < Nx/2)
            kx[i]=2*M_PI*(double)(i)/(double)(Nx*dx);
        else
            kx[i]=2*M_PI*(double)(i-Nx)/(double)(Nx*dx);
    }

    for(int j=0; j<Ny; j++){
        if(j<Ny/2)
            ky[j]=2*M_PI*(double)(j)/(double)(Ny*dy);
        else
            ky[j]=2*M_PI*(double)(j-Ny)/(double)(Ny*dy);
    }

    int ii;
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            ii=i*Nx+j;
            k2[ii]=(kx[i]*kx[i])+(ky[j]*ky[j]);
            k4[ii]=k2[ii]*k2[ii];
        }
    }
}
