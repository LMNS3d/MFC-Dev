# left files in pre_process/

## m_check_ib_patches.fpp

check user-provided inputs when using IB method, not utilized by our simulation.

## m_check_patches.fpp

check user-inputed parameters in each patch for different geometry, whether they supporting smoothing or not.

## m_checker.f90

check whether the inputs and the combination of the inputs are supported or not.

## m_global_parameters.fpp

define all the globally utilized variables and parameters.

### s_assign_default_values_to_user_inputs

...

### s_initialize_global_parameters_module

- model_eqns = 2 (five-equation model)
  - $\alpha_1\rho_1$: cont_idx%beg = 1
  - $\alpha_2\rho_2$: cont_idx%end = num_fluids
  - $\rho u$: mom_idx%beg = cont_idx%end + 1
  - $\rho v$: mom_idx%end = cont_idx%end + num_dims
  - $E$: E_idx = mom_idx%end + 1
  - $\alpha_1$: adv_idx%beg = E_idx + 1
  - $\alpha_2$: adv_idx%end = E_idx + num_fluids

sys_size = adv_idx%end

- model_eqns = 3 (six-equation model)
  - $\alpha_1\rho_1$: cont_idx%beg = 1
  - $\alpha_2\rho_2$: cont_idx%end = num_fluids
  - $\rho u$: mom_idx%beg = cont_idx%end + 1
  - $\rho v$: mom_idx%end = cont_idx%end + num_dims
  - $E$: E_idx = mom_idx%end + 1
  - $\alpha_1$: adv_idx%beg = E_idx + 1
  - $\alpha_2$: adv_idx%end = E_idx + num_fluids
  - $\alpha_1\rho_1 e_1$: internalEnergies_idx%beg = adv_idx%end + 1
  - $\alpha_2\rho_2 e_2$: internalEnergies_idx%end = adv_idx%end + num_fluids

sys_size = internalEnergies_idx%end

## m_grid.f90

create the rectilinear grid, which may be either uniform or non-uniform, see Johsen (2007) for non-uniform, which use hyperbolic tangent function.

### s_initialize_grid_module

point `s_generate_grid` to

- `s_generate_parallel_grid` for parallel_io = T
- `s_generate_serial_grid` otherwise

### s_finalize_grid_module

point `s_generate_grid` => null()

### s_generate_serial_grid

$$
dx = \frac{x_{end} - x_{beg}}{m+1} \\
x_{cc}(i) = x_{beg} + \frac{2i+1}{2}dx \\
x_{cb}(i-1) = x_{beg} + idx
$$

stretch the grid if needed.

### s_generate_parallel_grid

MPI grid generation

...

## m_model.fpp

read a binary STL file imported outside for mesh generation.

## m_mpi_proxy.fpp

Only processor with rank 0 is in charge of reading and checking user provided inputs, then broadcasting the required information.

### s_mpi_bcast_user_inputs

broadcast information from processor 0 to other ones.

### s_mpi_decompose_computational_domain

decompose computational domain for different processors.

## p_main.f90

1. read grid
2. apply initial condition
3. save data
