//
// Created by Jitin Nair on 14/01/20.
//

#include "spinodal_elasticity.h"


double* green_tensor(Nx, Ny, kx, ky, cm11, cm12, cm44, cp11, cp12, cp44, tmatrix){
    double c11, c12, c44, chi;
    c11 = 0.5*(cm11+cp11);
    c12 = 0.5*(cm12+cp12);
    c44 = 0.5*(cm44+cm44);
    chi=(c11-c12-2.0*c44)/c44;

    
    return tmatrix;
}