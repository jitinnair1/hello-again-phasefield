#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void write_to_VTK( int nx, int ny, int nz, double dx, double dy, double dz, int iprint, double conc[] );
