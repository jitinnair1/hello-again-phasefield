#include "set_array_zero.h"

/*
Output:

100
[1]    68713 segmentation fault  ./a.out

*/

int main(){
  int Nx=10, Ny=10;
  float array1[Nx][Ny];

  // arraysize of array1 before
  long arraySize1 = sizeof(array1)/sizeof(array1[0][0]);
  printf("%ld\n", arraySize1);

  // Pass to function to initialize
  set_array_zero(Nx, Ny, array1);

  // arraysize of array1 after
  long arraySize2 = sizeof(array1)/sizeof(array1[0][0]);
  printf("%ld\n", arraySize2);

  float array_new[Nx][Ny];
  set_array_zero(Nx, Ny, array_new);

}
