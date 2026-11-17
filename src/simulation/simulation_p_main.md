# left files in simulation

## p_main.fpp

time-stepping loop

1. `s_perform_time_step` to do simulation for one time step

2. `s_save_data` to save data when
    - (t_step - t_step_start, t_step_save) = 0
    - t_step = t_step_stop

3. `s_save_performance_metrics` to save time information at t_step_stop

## m_derived_variables.f90

to derive the unadvected volume fraction, specific heat ratio, liquid stiffness,
speed of sound, vorticity and the numerical Schlieren function

But the calculation process of these variables not found in this file.

### s_initialize_derived_variables

computing FD coefficients. `s_compute_finite_difference_coefficients`

### s_compute_derived_variables

only use `s_derive_acceleration_component` with RK sub-step solutions to calculate
acceleration component

## m_checker.fpp

check for compatible input files.

## m_global_parameters.fpp

Global Variables:

- proc_rank: Rank of the local processor
- x_cb: cell-boundary location
- x_cc: cell-center location
- dx
- buff_size:
  - number of cells necessary to store enough boundary conditions data
- start_idx: starting cell-center index of local processor in global grid
  - for parallel computation, grid separation.
- gammas: $\Gamma_k$
- gs_min: $\gamma_k = \frac{1}{\Gamma_k} + 1$
- pi_infs: $\Pi_k$
- ps_inf: $p^{\infty}_k = \frac{\Pi_k}{\Gamma_k + 1}$
- cvs: $C_{v,k}$
- qvs: $q_k$
- qvps: $q_k'$

### s_assign_default_values_to_user_inputs

Assign default values to the user inputs before reading them in.

### s_initialize_global_parameters_module

$$
weno_{polyn} = \frac{weno_{order} - 1}{2}
$$

- model_eqns = 1 (Gamma/Pi_inf Model)

  - num_fluids = 1
  - cont_idx%beg = 1
  - cont_idx%end = cont_idx%beg
  - mom_idx%beg = cont_idx%end + 1
  - mom_idx%end = cont_idx%end + num_dims
  - E_idx = mom_idx%end + 1
  - adv_idx%beg = E_idx + 1
  - adv_idx%end = adv_idx%beg + 1
  - gamma_idx = adv_idx%beg
  - pi_inf_idx = adv_idx%end

  sys_size = adv_idx%end

- model_eqns = 2 (five-equation mdoel)

  - cont_idx%beg = 1
  - cont_idx%end = num_fluids
  - mom_idx%beg = cont_idx%end + 1
  - mom_idx%end = cont_idx%end + num_dims
  - E_idx = mom_idx%end + 1
  - adv_idx%beg = E_idx + 1
  - adv_idx%end = E_idx + num_fluids

  sys_size = adv_idx%end

- model_eqns = 3 (six-equation model)

  - cont_idx%beg = 1
  - cont_idx%end = num_fluids
  - mom_idx%beg = cont_idx%end + 1
  - mom_idx%end = cont_idx%end + num_dims
  - E_idx = mom_idx%end + 1
  - adv_idx%beg = E_idx + 1
  - adv_idx%end = E_idx + num_fluids
  - internalEnergies_idx%beg = adv_idx%end + 1
  - internalEnergies_idx%end = adv_idx%end + num_fluids

  sys_size = internalEnergies_idx%end

$$
buff_{size} = weno_{polyn} + 2 \\
startx = -buff_{size} \\
n>0: starty = -buff_{size} \\
p>0: startz = -buff_{size}
$$

- momxb = mom_idx%beg
- momxe = mom_idx%end
- advxb = adv_idx%beg
- advxe = adv_idx%end
- contxb = cont_idx%beg
- contxe = cont_idx%end
- intxb = internalEnergies_idx%beg
- intxe = internalEnergies_idx%end
- x_cb, y_cb (n>0), z_cb (p>0)
  - -1-buff_size : m/n/p+buff_size
- x_cc, y_cc (n>0), z_cc (p>0)
  - -buff_size : m/n/p+buff_size
- dx, dy (n>0), dz (p>0)
  - -buff_size : m/n/p+buff_size

### s_finalize_global_parameters_module

deallocate x_cc, x_cb, dx, start_idx.
