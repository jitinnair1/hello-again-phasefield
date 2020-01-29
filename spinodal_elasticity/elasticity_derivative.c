//
// Created by Jitin Nair on 24/01/20.
//

#include "spinodal_elasticity.h"


void elasticity_derivative(int Nx, int Ny, double tmatx, double kx, double ky,
                           double s11, double s22, double s12,double e11, double e22, double e12,
                           double ed11 ,double ed22 ,double ed12 ,double cm11 ,double cm12 ,double cm44 ,
                           double cp11 ,double cp12 ,double cp44 ,double ea, double ei0, fftw_complex* conc,
                           fftw_plan p4, fftw_plan p10 ) {

    int niter=10,
    num_points = Nx * Ny;

    double tolerance=0.001;
    double sum_norm=0.0,
    old_norm=0.0,
    conver=0.0,
    normF;

    for (int ii = 0; ii < num_points; ++ii) {
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
            for (int j = 0; j < Ny; ++j) {
                smatx[i][j][1][1] = s11k;
                smatx[i][j][1][2] = s12k;
                smatx[i][j][2][1] = s12k;
                smatx[i][j][2][2] = s22k;

                ematx[i][j][1][1] = e11k;
                ematx[i][j][1][2] = e12k;
                ematx[i][j][2][1] = e12k;
                ematx[i][j][2][2] = e22k;

                for (int ii = 0; ii < 2; ++ii) {
                    for (int jj = 0; jj < 2; ++jj) {
                        for (int kk = 0; kk < 2; ++kk) {
                            for (int ll = 0; ll < 2; ++ll) {
                                ematx[i][j][ii][jj] =
                                        ematx[i][j][ii][jj] - tmatx[i][j][ii][jj][kk][ll] * smatx[i][j][ii][jj];
                            }
                        }
                    }
                }

                //get updated values of strains in fourier space
                e11k = ematx[i][j][1][1];
                e22k = ematx[i][j][2][2];
                e12k = ematx[i][j][1][2];
            }
        }

        //take strains to real space
        fftw_execute_dft(p10, e11k, e11);
        fftw_execute_dft(p10, e12k, e12);
        fftw_execute_dft(p10, e22k, e22);

        //calculate stress in real space
        sum_norm = 0.0;
        for (int ii = 0; ii < num_points; ++ii) {
            s11[ii]=c11[ii]*(ea[0]+e11[ii]-ei11[ii]-ed11[ii])+c12[ii]*(ea[1]+e22[ii]-ei22[ii]-ed22[ii]);
            s22[ii]=c12[ii]*(ea[1]+e22[ii]-ei22[ii]-ed22[ii])+c12[ii]*(ea[0]+e11[ii]-ei11[ii]-ed11[ii]);
            s12[ii]=2.0*c22[ii]*(ea[2]+e12[ii]-ei12[ii]-ed12[ii]);
            sum_stres[ii] = s11[ii]+s22[ii]+s12[ii];
            sum_norm=sum_norm+(sum_stres[ii]*sum_stres[ii]);
        }

        //get euclidean norm
        normF = sqrt(sum_norm);

        //check for convergence
        if(k != 1)
            conver = fabs((normF-old_norm)/old_norm);
        if(conver <= tolerance)
            break;
        old_norm = normF;
    }

    //return value of derivative
    for (int ii = 0; ii < num_points; ++ii) {
        et11[ii] =ea[0]+e11[ii]-ei11[ii]-ed11[ii];
        et22[ii] =ea[1]+e22[ii]-ei22[ii]-ed22[ii];
        et12[ii] =ea[2]+e12[ii]-ei12[ii]-ed12[ii];


        delsdc[ii] = 0.5*(et11[ii]*((cp12-cm12)*et22[ii]+(cp11-cm11)*et11[ii]-c12[ii]*ei0-c11[ii]*ei0)-ei0*(c12[ii]*et22[ii]
                 +c11[ii]*et11[ii]) +  ((cp11-cm11)*et22[ii]+(cp12-cm12)*et11[ii]-c12[ii]*ei0-c11[ii]*ei0)*et22[ii]
        -ei0*(c11[ii]*et22[ii]+c12[ii]*et11[ii]) + 2.0*(cp44-cm44)*et12[ii]*et12[ii]-4.0*ei0*c44[ii]*et12[ii]);
    }

}



