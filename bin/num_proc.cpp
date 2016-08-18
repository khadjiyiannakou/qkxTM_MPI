#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(){

  printf("You use a function with omp support\n");
  int n_cores = omp_get_num_procs();
  printf("The system you are running has %d cores\n",n_cores);
  int max_n_threads = omp_get_max_threads();
  printf("You asked for %d threads\n",max_n_threads);

  return 0;
}
