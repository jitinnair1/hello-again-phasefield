#include "spinodal_spect.h"

float** array_allocate(int nx, int ny, float **array){
  int m=nx;

  // dynamically create array of pointers of size m
  array = calloc(m, sizeof(float *));

  // dynamically allocate memory of size n for each row
  for (int r = 0; r <= ny; r++)
  array[r] = calloc(ny, sizeof(float));

  return array;
}
