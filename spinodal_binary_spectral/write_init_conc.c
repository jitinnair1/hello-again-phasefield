//
// Created by Jitin Nair on 13/01/20.
//
#include "spinodal_spect.h"

void write_init_conc(int Nx, int Ny, int NxNy, double conc_print[NxNy]){
    FILE *file;
    char filename[30];
    sprintf(filename, "./output/conc0FFTWbinary.txt");
    file = fopen(filename,"w");
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; ++j) {
            fprintf(file, "%f", conc_print[i*Ny+j]);
            fprintf(file,"\t");
        }
        fprintf(file,"\n");
    }
    fclose(file);
}
