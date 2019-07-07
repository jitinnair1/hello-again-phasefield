#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void write_to_VTK( int nx, int ny, int nz, float dx, float dy, float dz, int iprint, float conc[] );
