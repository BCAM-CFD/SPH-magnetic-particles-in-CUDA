/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

#include "class_system.h"
#include "config.h"
#include "kernel_functions.h"
#include <stdio.h>

// Function to do the second step of the Velocity Verlet algorithm
int class_system::V_Verlet_step2(dim3 numBlocks,
				 dim3 threadsPerBlock,
				 real* k_vx,
				 real* k_vy,
				 real* k_vz,
				 real* k_fx,
				 real* k_fy,
				 real* k_fz,
				 real* k_mass,
				 int*  k_type,
				 real* k_coll_vx,
				 real* k_coll_vy,
				 real* k_coll_vz,
				 real* k_coll_omegax,
				 real* k_coll_omegay,
				 real* k_coll_omegaz,				 
				 real* k_fx_colloid,
				 real* k_fy_colloid,
				 real* k_fz_colloid,
				 real* k_tx_colloid,
				 real* k_ty_colloid,
				 real* k_tz_colloid,				 
				 int   coll_move) {

  cudaError_t cuda_err;

  //--- Fluid and walls particles are moved ---      
  kernel_move_sys_VV_part2<<<numBlocks, threadsPerBlock>>>(k_vx, k_vy, k_vz,
							   k_fx, k_fy, k_fz, k_mass, k_type);
  cuda_err = cudaGetLastError();
  if (cuda_err != cudaSuccess) {
    printf("Error in kernel_move_sys: %s\n", cudaGetErrorString(cuda_err));
    return 1;
  }

  //--- Particles of colloids are moved ---
  if (coll_move == 0) 
    if (N_colloids > 0) {
      kernel_move_colloids_VV_part2<<<numBlocks, threadsPerBlock>>>(k_coll_vx,
								    k_coll_vy,
								    k_coll_vz,
								    k_coll_omegax,
								    k_coll_omegay,
								    k_coll_omegaz,
								    k_fx_colloid,
								    k_fy_colloid,
								    k_fz_colloid,
								    k_tx_colloid,
								    k_ty_colloid,
								    k_tz_colloid);
      cuda_err = cudaGetLastError();
      if (cuda_err != cudaSuccess) {
	printf("Error in kernel_colloids_VV_part1: %s\n", cudaGetErrorString(cuda_err));
	return 1;
      }
    }
  
  cudaDeviceSynchronize();  // We require the kernel to end to continue    

  return 0;

}
