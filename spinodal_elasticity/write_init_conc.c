//
// Created by Jitin Nair on 13/01/20.
//
#include "spinodal_elasticity.h"

void write_init_conc(int Nx, int Ny, int NxNy, double conc_print[NxNy], char label[30]){

    FILE *file;
    char fileToOpen[strlen(label) + 5];

    strcpy(fileToOpen, label);
    strcat(fileToOpen, ".txt");

    file = fopen(fileToOpen,"w");

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; ++j) {
            fprintf(file, "%f", conc_print[i*Ny+j]);
            fprintf(file,"\t");
        }
        fprintf(file,"\n");
    }
    fclose(file);
}
