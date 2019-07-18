#include "set_array_zero.h"

/*
Output:

Array size before is : 100
Array size after is: 0

*/

int main(){
  int Nx=10, Ny=10;
  float array1[Nx][Ny];

  // arraysize of array1 before
  long arraySize1 = sizeof(array1)/sizeof(array1[0][0]);
  printf("Array size before is : %ld\n", arraySize1);

  // Pass to function to initialize
  set_array_zero(Nx, Ny, array1);

  // arraysize of array1 after
  long arraySize2 = sizeof(array1)/sizeof(array1[0][0]);
  printf("Array size after is: %ld\n", arraySize2);

}
