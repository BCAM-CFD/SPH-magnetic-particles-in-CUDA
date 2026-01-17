/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "kernel_functions.h"
#include "config.h"

//Part 2 of the function to move particles with a velocity Verlet method with lambda = 0.5.
__global__ void kernel_move_sys_VV_part2(real* __restrict__ vx,
					 real* __restrict__ vy,
					 real* __restrict__ vz,				
					 real* __restrict__ fx,
					 real* __restrict__ fy,
					 real* __restrict__ fz,
					 real* __restrict__ mass,
					 int* __restrict__ type) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= N) return;

  if (type[i] == 0) { //  If it is from fluid
    real half_dt_over_mass = 0.5 * dt / mass[i];    
    // Velocity at t + dt/2 
    vx[i] = vx[i] + half_dt_over_mass * fx[i];
    vy[i] = vy[i] + half_dt_over_mass * fy[i];
    if (dim == 3)
      vz[i] = vz[i] + half_dt_over_mass * fz[i];
    
  }
}

