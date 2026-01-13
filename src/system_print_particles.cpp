#include <iostream>
#include "cuda_runtime.h"
#include "class_system.h"

// function to print system info
void class_system::print_particles(real* k_x,
				   real* k_y,
				   real* k_z,
				   real* k_vx,
				   real* k_vy,
				   real* k_vz,
				   real* k_mass,				   
				   real* k_dens,
				   real* k_press,				   
				   int step) {
  int Npart = this->N;

  // Note that type does not change, so we do not need to copy its value from the kernel to the host

  char filename[50];
  sprintf(filename, "micro-%d.dat", step);
  FILE* file = fopen(filename, "w");
  if (!file) {
    printf("System print particle error: Error opening the file %s\n", filename);
    return;
  }
  
  cudaMemcpy(this->x,     k_x,     Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->y,     k_y,     Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->z,     k_z,     Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->vx,    k_vx,    Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->vy,    k_vy,    Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->vz,    k_vz,    Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->mass,  k_mass,  Npart * sizeof(real), cudaMemcpyDeviceToHost);  
  cudaMemcpy(this->dens,  k_dens,  Npart * sizeof(real), cudaMemcpyDeviceToHost);
  cudaMemcpy(this->press, k_press, Npart * sizeof(real), cudaMemcpyDeviceToHost);
  if (dim == 2)
    for (int i = 0; i < Npart; ++i)
      fprintf(file, "%d " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " %d\n",
	      i,
	      this->x[i],
	      this->y[i],
	      this->vx[i],
	      this->vy[i],
	      this->mass[i],	      
	      this->dens[i],
	      this->press[i],
	      this->type[i]);
  else //--- dim == 3  ---
    for (int i = 0; i < Npart; ++i)
      fprintf(file, "%d " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " " REAL_FMT " %d\n",
	      i,
	      this->x[i],
	      this->y[i],
	      this->z[i],
	      this->vx[i],
	      this->vy[i],
	      this->vz[i],
	      this->mass[i],	      
	      this->dens[i],
	      this->press[i],
	      this->type[i]);    
  fclose(file);

  printf("Micro file written. Time step %d\n", step);
}
