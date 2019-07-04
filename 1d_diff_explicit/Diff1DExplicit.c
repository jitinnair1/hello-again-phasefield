#include "stdio.h"
#include "math.h"

int main() {

  //Declaring variables
  float alpha, D, dx, dt;
  int j, k, i, istep, N;
  FILE *fpointer;
  fpointer=fopen("profile.txt", "w");

  //Initial system configuration
  N=128;
  D=1.0;
  dx=0.5;
  dt=0.1;
  istep=1500;
  alpha=D*dt/(dx*dx);

  //Initial conditions and profile
  float conc[N];
  conc[0]=1;
    for (i=1; i<=N; i++){
        conc[i]=0;
    }

  //Evolution of profile

  //loop for saving results
  for (k=1; k<8; k++) {

    //loop for evolution
    for (j=1; j<istep; j++) {

      //loop for explicit condition
      for (i = 1; i < N; i++) {
        conc[i]=(1-2*alpha)*conc[i]+alpha*(conc[i-1]+conc[i+1]);
      }

      //handling boundary condition at right end
      conc[N]=(1-2*alpha)*conc[N]+2*alpha*(conc[N-1]);
    }

    //Write results to file in fpointer after istep iterations
    for (i = 0; i < N; i++) {
      fprintf(fpointer, "%f\n", conc[i]);
    }

    fprintf(fpointer, "\n");
  }
  fclose(fpointer);

}
