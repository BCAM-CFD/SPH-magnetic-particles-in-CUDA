#include "class_system.h"
#include "cuda_runtime.h"

#include <stdio.h>

// Function to copy the system data to the device
void class_system::copy_pointers_to_device(real** k_x,
					   real** k_y,
					   real** k_z,
					   real** k_vx,
					   real** k_vy,
					   real** k_vz,
					   real** k_mass,
					   real** k_fx,
					   real** k_fy,
					   real** k_fz,
					   real** k_dens,
					   real** k_press,
					   int**  k_particle_cell,
					   int**  k_particle_index,
					   int**  k_cell_start,
					   int**  k_cell_end,
					   real** k_kin_energy,
					   int**  k_type,
					   real** k_fx_wall,
					   real** k_fy_wall,
					   real** k_fz_wall,
					   int**  k_walls_list,
					   int**  k_walls_start,
					   int**  k_colloids_list,
					   int**  k_colloids_start,
					   real** k_coll_x,
					   real** k_coll_y,
					   real** k_coll_z,
					   real** k_coll_vx,
					   real** k_coll_vy,
					   real** k_coll_vz,
					   real** k_coll_omegax,
					   real** k_coll_omegay,
					   real** k_coll_omegaz,
					   real** k_coll_theta,
					   real** k_coll_quat0,
					   real** k_coll_quatx,
					   real** k_coll_quaty,
					   real** k_coll_quatz,
					   real** k_x_center,
					   real** k_y_center,
					   real** k_z_center,
					   real** k_fx_colloid,
					   real** k_fy_colloid,
					   real** k_fz_colloid,
					   real** k_tx_colloid,
					   real** k_ty_colloid,
					   real** k_tz_colloid,
					   int**  k_coll_cell,
					   int**  k_coll_index,
					   int**  k_coll_cell_start,
					   int**  k_coll_cell_end,
					   real** k_magnetic_mom) {
  
  //------- Variables to copy -------------
  cudaMalloc(k_x             , N * sizeof(real));
  cudaMalloc(k_y             , N * sizeof(real));
  cudaMalloc(k_z             , N * sizeof(real));
  cudaMalloc(k_vx            , N * sizeof(real));
  cudaMalloc(k_vy            , N * sizeof(real));
  cudaMalloc(k_vz            , N * sizeof(real));
  cudaMalloc(k_mass          , N * sizeof(real));
  cudaMalloc(k_kin_energy    , sizeof(real));
  cudaMalloc(k_type          , N * sizeof(int));
  if (wall) {
    cudaMalloc(k_walls_list    , Nlist_walls * sizeof(int));
    cudaMalloc(k_walls_start   , Nwalls * sizeof(int));
  }
  cudaMalloc(k_colloids_list , Nlist_colloids * sizeof(int));
  cudaMalloc(k_colloids_start, N_colloids * sizeof(int));  
  cudaMalloc(k_coll_x        , N_colloids * sizeof(real));
  cudaMalloc(k_coll_y        , N_colloids * sizeof(real));
  cudaMalloc(k_coll_z        , N_colloids * sizeof(real));
  cudaMalloc(k_coll_vx       , N_colloids * sizeof(real));
  cudaMalloc(k_coll_vy       , N_colloids * sizeof(real));
  cudaMalloc(k_coll_vz       , N_colloids * sizeof(real));
  cudaMalloc(k_coll_omegax   , N_colloids * sizeof(real));
  cudaMalloc(k_coll_omegay   , N_colloids * sizeof(real));
  cudaMalloc(k_coll_omegaz   , N_colloids * sizeof(real));
  cudaMalloc(k_coll_theta    , N_colloids * sizeof(real));
  cudaMalloc(k_coll_quat0    , N_colloids * sizeof(real));
  cudaMalloc(k_coll_quatx    , N_colloids * sizeof(real));
  cudaMalloc(k_coll_quaty    , N_colloids * sizeof(real));
  cudaMalloc(k_coll_quatz    , N_colloids * sizeof(real));        
  cudaMalloc(k_x_center      , Nlist_colloids * sizeof(real));
  cudaMalloc(k_y_center      , Nlist_colloids * sizeof(real));
  cudaMalloc(k_z_center      , Nlist_colloids * sizeof(real));

  cudaMemcpy(*k_x,     x,     N * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_y,     y,     N * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_z,     z,     N * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_vx,    vx,    N * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_vy,    vy,    N * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_vz,    vz,    N * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_mass,  mass,  N * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_kin_energy, &kin_energy,  sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_type,  type,  N * sizeof(int), cudaMemcpyHostToDevice);
  if (wall) {
    cudaMemcpy(*k_walls_list, walls_list, Nlist_walls * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(*k_walls_start, walls_start, Nwalls * sizeof(int), cudaMemcpyHostToDevice);
  }
  cudaMemcpy(*k_colloids_list, colloids_list,
	     Nlist_colloids * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_colloids_start, colloids_start,
	     N_colloids * sizeof(int), cudaMemcpyHostToDevice);  
  cudaMemcpy(*k_coll_x,      coll_x     , N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_coll_y,      coll_y     , N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_coll_z,      coll_z     , N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_coll_vx,     coll_vx    , N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_coll_vy,     coll_vy    , N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_coll_vz,     coll_vz    , N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_coll_omegax, coll_omegax, N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_coll_omegay, coll_omegay, N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_coll_omegaz, coll_omegaz, N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_coll_theta,  coll_theta , N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_coll_quat0,  coll_quat0 , N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_coll_quatx,  coll_quatx , N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_coll_quaty,  coll_quaty , N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_coll_quatz,  coll_quatz , N_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_x_center,   x_center    , Nlist_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_y_center,   y_center    , Nlist_colloids * sizeof(real), cudaMemcpyHostToDevice);
  cudaMemcpy(*k_z_center,   z_center    , Nlist_colloids * sizeof(real), cudaMemcpyHostToDevice);

  //--------- Other variables --------------
  cudaMalloc(k_fx             , N * sizeof(real));
  cudaMalloc(k_fy             , N * sizeof(real));
  cudaMalloc(k_fz             , N * sizeof(real));
  cudaMalloc(k_dens           , N * sizeof(real));
  cudaMalloc(k_press          , N * sizeof(real));    
  cudaMalloc(k_particle_cell  , N * sizeof(int));
  cudaMalloc(k_particle_index , N * sizeof(int));
  cudaMalloc(k_cell_start     , Ntotal_cells * sizeof(int));
  cudaMalloc(k_cell_end       , Ntotal_cells * sizeof(int));
  if (wall) {
    cudaMalloc(k_fx_wall        , Nwalls * sizeof(real));
    cudaMalloc(k_fy_wall        , Nwalls * sizeof(real));
    cudaMalloc(k_fz_wall        , Nwalls * sizeof(real));
  }
  cudaMalloc(k_fx_colloid     , N_colloids * sizeof(real));
  cudaMalloc(k_fy_colloid     , N_colloids * sizeof(real));
  cudaMalloc(k_fz_colloid     , N_colloids * sizeof(real));
  cudaMalloc(k_tx_colloid     , N_colloids * sizeof(real));
  cudaMalloc(k_ty_colloid     , N_colloids * sizeof(real));
  cudaMalloc(k_tz_colloid     , N_colloids * sizeof(real));
  cudaMalloc(k_coll_cell      , N_colloids * sizeof(int));
  cudaMalloc(k_coll_index     , N_colloids * sizeof(int));
  cudaMalloc(k_coll_cell_start, Ntotal_cells_colloids * sizeof(int));
  cudaMalloc(k_coll_cell_end  , Ntotal_cells_colloids * sizeof(int));
  cudaMalloc(k_magnetic_mom   , 3 * sizeof(real));  
    
  cudaMemset(*k_fx             , 0, N * sizeof(real));
  cudaMemset(*k_fy             , 0, N * sizeof(real));
  cudaMemset(*k_fz             , 0, N * sizeof(real));
  cudaMemset(*k_dens           , 0, N * sizeof(real));              
  cudaMemset(*k_press          , 0, N * sizeof(real));
  cudaMemset(*k_particle_cell  , 0, N * sizeof(int));              
  cudaMemset(*k_particle_index , 0, N * sizeof(int));
  cudaMemset(*k_cell_start     , 0, Ntotal_cells * sizeof(int));              
  cudaMemset(*k_cell_end       , 0, Ntotal_cells * sizeof(int));
  if (wall) {
    cudaMemset(*k_fx_wall        , 0, Nwalls * sizeof(real));
    cudaMemset(*k_fy_wall        , 0, Nwalls * sizeof(real));
    cudaMemset(*k_fz_wall        , 0, Nwalls * sizeof(real));
  }
  if (N_colloids > 0) {
    cudaMemset(*k_fx_colloid     , 0, N_colloids * sizeof(real));
    cudaMemset(*k_fy_colloid     , 0, N_colloids * sizeof(real));
    cudaMemset(*k_fz_colloid     , 0, N_colloids * sizeof(real));
    cudaMemset(*k_tx_colloid     , 0, N_colloids * sizeof(real));
    cudaMemset(*k_ty_colloid     , 0, N_colloids * sizeof(real));
    cudaMemset(*k_tz_colloid     , 0, N_colloids * sizeof(real));
    cudaMemset(*k_coll_cell      , 0, N_colloids * sizeof(int));              
    cudaMemset(*k_coll_index     , 0, N_colloids * sizeof(int));
    cudaMemset(*k_coll_cell_start, 0, Ntotal_cells_colloids * sizeof(int));              
    cudaMemset(*k_coll_cell_end  , 0, Ntotal_cells_colloids * sizeof(int));
    cudaMemset(*k_magnetic_mom   , 0, 3 * sizeof(real));
  }

}
