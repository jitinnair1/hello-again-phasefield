#include "set_array_zero.h"

/*
Output:

For array1:
Array size before is : 100
Array size after is: 100
For array2:
Array size before is : 100
Array size after is : 100

*/

int main(){
  int Nx=10, Ny=10;
  float array1[Nx][Ny];

  printf("For array1:\n");

  // arraysize of array1 before
  long arraySize1 = sizeof(array1)/sizeof(array1[0][0]);
  printf("Array size before is : %ld\n", arraySize1);

  // Pass to function to initialize
  set_array_zero(Nx, Ny, array1);

  // arraysize of array1 after
  long arraySize2 = sizeof(array1)/sizeof(array1[0][0]);
  printf("Array size after is: %ld\n", arraySize2);

  printf("For array2:\n");
  float array2[Nx][Ny];

  // arraysize of array2 before
  long arraySize3 = sizeof(array2)/sizeof(array2[0][0]);
  printf("Array size before is : %ld\n", arraySize3);

  // Pass to function to initialize
  set_array_zero(Nx, Ny, array2);

  // arraysize of array2 before
  long arraySize4 = sizeof(array2)/sizeof(array2[0][0]);
  printf("Array size after is : %ld\n", arraySize4);

}
