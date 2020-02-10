// Created by Jitin Nair on 14/01/20.

#include "spinodal_elasticity.h"

double *green_tensor(int Nx, int Ny, double kx[], double ky[], double cm11, double cm12,
                     double cm44, double cp11, double cp12, double cp44, double tmatx[][Ny][2][2][2][2]){

    double c11, c12, c44, chi;
    double rr, d0;
    double omeg11[Nx][Ny], omeg22[Nx][Ny], omeg12[Nx][Ny];
    double gmatx[2][2];
    double dvect[2];

    c11 = 0.5*(cm11+cp11);
    c12 = 0.5*(cm12+cp12);
    c44 = 0.5*(cm44+cm44);
    chi=(c11-c12-2.0*c44)/c44;

    //Get omeg_mat values
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {

            // calculate d0
            rr=kx[i]*kx[i]+ky[j]*ky[j];
            d0=c11*rr*rr*rr+chi*(c11+c12)*rr*(kx[i]*kx[i]*ky[j]*ky[j]);

            //check for small rr
            if(rr < 1.0e-8)
                d0=1.0;

            //calculate omeg_mat values
            omeg11[i][j]=(c44*rr*rr+(c11-c44)*rr*ky[j]*ky[j])/(c44*d0);
            omeg22[i][j]=(c44*rr*rr+(c11-c44)*rr*kx[i]*kx[i])/(c44*d0);
            omeg12[i][j]=-(c12+c44)*kx[i]*ky[j]*rr/(c44*d0);

        }
    }

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            //green tensor
            gmatx[0][0]=omeg11[i][j];
            gmatx[0][1]=omeg12[i][j];
            gmatx[1][0]=omeg12[i][j];
            gmatx[1][1]=omeg22[i][j];

            //position vector
            dvect[0]=kx[i];
            dvect[1]=ky[j];

            //green operator
            for (int ii = 0; ii < 2; ++ii) {
                for (int jj = 0; jj < 2; ++jj) {
                    for (int kk = 0; kk < 2; ++kk) {
                        for (int ll = 0; ll < 2; ++ll) {
                            tmatx[i][j][ii][jj][kk][ll] = 0.25 * (
                                    gmatx[jj][kk] * dvect[ll] * dvect[ii]
                                    + gmatx[ii][kk] * dvect[ll] * dvect[jj]
                                    + gmatx[jj][ll] * dvect[kk] * dvect[ii]
                                    + gmatx[ii][ll] * dvect[kk] * dvect[jj]);
                        }
                    }
                }
            }
        }
    }
    
    return tmatx;
}
