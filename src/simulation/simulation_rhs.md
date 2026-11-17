# m_rhs.fpp

This code is generally one of the three most important files, while the other two
are m_weno.fpp, m_riemann_solver.fpp, m_rhs.fpp.

Variables:

- q_cons_qp:
  - vector_field
  - conservative variables at cell-interior Gaussian quadrature points (QP),
  - recosntructed from the cell-average conservative variables q_cons_vf with WENO scheme
- q_prim_qp:
  - vector_field
  - primitive variables at cell-interior QP,
  - calculated from q_cons_qp
- flux_n:
  - (:) of vector_field
  - cell-boundary values of the fluxes
  - use Riemann solver based on the left and right cell-boundary values of the primitive variables
- flux_src_n:
  - (:) of vector_field
  - source term
  - based on Riemann solver
- qL_prim, qR_prim:
  - (:) of vector field
  - left and right cell-boundary reconstructed variables
  - use `s_weno` to calculate
- qL_rsx_vf, qR_rsx_vf, qL_rxy_vf, qR_rsy_vf, ...
  - (:,:,:,:)
  - reshaped for dimension-splitting
- blkmod1, blkmod2, alpha1, alpha2, Kterm
  - (:,:,:)
  - variables for calculation of K term in Kapila's five-equation model
- gamma_min, pres_inf
  - (:)

## s_initialize_rhs_module

memory allocation, parameters computation.

- Coordinate direction indexes
  - ix%beg = -buff_size
  - ix%end = m - ix%beg
  - similar for other directions
- allocate q_cons_qp, q_prim_qp
  - for $\alpha_k\rho_k,\alpha_k$, share the same space
  - only one point within one cell
- allocate qL_prim, qR_prim
  - only allocate spatial index for mom_idx
- allocate qL_rsx_vf, qR_rsx_vf, ...
  - `(:,:,:,:)`: the first three are spatial indexes, the last one is equation index
  - reshaped for different directions
- allocate flux_n, flux_src_n
- calculate $\gamma_k,p^{\infty}_k$
  - $\gamma_k = \frac{1}{\Gamma_k} + 1$
  - $p^{\infty}_k = \frac{\Pi_k}{\Gamma_k + 1}$
- point s_riemann_solver to
  - s_hll_riemann_solver (1)
  - s_hllc_riemann_solver (2)
- point s_convert_to_mixture_variables to
  - s_convert_mixture_to_mixture_variables (gamma/pi_inf model)
  - s_convert_species_to_mixture_variables (volume fraction model)

## s_compute_rhs

Variables:

- q_cons_vf:
  - (sys_size) of scalar_field
  - cell-average conservative variables
  - inout
- q_prim_vf:
  - (sys_size) of scalar_field
  - cell-average primitive variables
  - inout
- rhs_vf:
  - (sys_size) of scalar_field
  - inout

1. Coordinate direction indexes

     - ix%beg = -buff_size
     - ix%end = m - ix%beg
     - similar for other directions

2. q_cons_qp = q_cons_vf

3. `s_convert_conservative_to_primitive_variables`
   - initially only the cell-average conservative variables are available
   - convert conservative to primitive cell-average variables to mitigate the oscillation around the interface

4. `s_populate_primitive_variables_buffers` to achieve the boundary values of the region according to different boundary conditions.

5. dimensional splitting loop

   1. coordinate direction indexes (for boundary value reconstruction)
     - ix%beg = -buff_size
     - ix%end = m - ix%beg
     - similar for other directions

   2. reconstruct cell-boundary primitive/conservative variables
     - `s_reconstruct_cell_boundary_values`
     - qL_rsx_vf, qR_rsx_vf
     - qL_rsy_vf, qR_rsy_vf
     - qL_rsz_vf, qR_rsz_vf

   3. coordinate direction indexes (for Flux calculation)
     - ix%beg = -1, iy%beg = 0, iz%beg = 0
     - ix%end = m, iy%end = n, iz%end = p
     - similar for other dimensions

   4. `s_riemann_solver`

   5. `s_compute_advection_source_term`

## s_compute_advection_source_term

Note that during this step, we do not see the parameter $q$, which is existed in the new stiffened-gas EOS, usually utilized in phase-change problem:

$$
\rho (e-q) = \Gamma p + \Pi
$$

Therefore, the part may need to be modified when considering cavitation problems.

Variables:

- q_cons_vf:
  - vector_field
  - cell-average conservative variables
  - inout
- q_prim_vf:
  - vector_field
  - cell-average primitive variables
  - inout
- rhs_vf:
  - (sys_size) of scalar_field
  - inout
- flux_src_n_vf:
  - vector_field
  - inout

1. calculate K term if the model is Kapila's five-equation model

$$
B_1 = blkmod_1 = \frac{(\Gamma_1 + 1)p + \Pi_1}{\Gamma_1} = \gamma_1 (p + p^{\infty}_1) \\
B_2 = blkmod_2 = \frac{(\Gamma_2 + 1)p + \Pi_2}{\Gamma_2} = \gamma_2 (p + p^{\infty}_2) \\
K = \frac{\alpha_1\alpha_2(B_2 - B_1)}{\alpha_1 B_2 + \alpha_2 B_1}
$$

2. calculate convective flux

- flux_n(k-1) corresponds to $F_{k-\frac{1}{2}}$
- flux_n(k) corresponds to $F_{k+\frac{1}{2}}$

$$
rhs = \frac{1}{dx}(F_{i-\frac{1}{2}} - F_{i+\frac{1}{2}})
$$

3. add source term for internal energy function of six-equation model ($\alpha_k\rho_k e_k$)

- flux_src_n(1)%vf(advxb)%sf(j, k, l) corresponds to $u_{j+\frac{1}{2}}$, which is calculated from HLLC Riemann solver
- flux_src_n(1)%vf(advxb)%sf(j-1, k, l) corresponds to $u_{j-\frac{1}{2}}$

$$
rhs = rhs - \frac{1}{dx}\alpha_kp(u_{i+\frac{1}{2}}-u_{i-\frac{1}{2}})
$$

4. add source term for advection equation

    - HLL

    how to calculate `flux_src_n`

    $$
    rhs = rhs + \frac{1}{dx}u (f_{src,i-\frac{1}{2}} - f_{src,i+\frac{1}{2}})
    $$

    - HLLC with K term (kapila's model)

    $$
    rhs_1 = rhs_1 + \frac{1}{dx}(\alpha_1 + K)(u_{i+\frac{1}{2}} - u_{i-\frac{1}{2}}) \\
    rhs_2 = rhs_2 + \frac{1}{dx}(\alpha_2 - K)(u_{i+\frac{1}{2}} - u_{i-\frac{1}{2}})
    $$

    - HLLC without K term (generally utilized)

    $$
    rhs = rhs + \frac{1}{dx}\alpha_k(u_{i+\frac{1}{2}} - u_{i-\frac{1}{2}})
    $$

## s_pressure_relaxation_procedure

utilized in six-equation model. After hyperbolic steps, the pressures of two phases are not equal to each other. This procedure is to infinitely relax the pressures from the internal-energy equations to a unique pressure.

Note that this procedure is not activated if phase change is activated.

Changed variables:

- volume fractions of each phase
- pressure

1. check whether need to relax or not
    - $\alpha_k > 1 - \epsilon$, if single phase, no need to relax

2. calculat pressure

    - initialization

    if $\alpha_k > \epsilon$

    $$
    p_k = (\frac{\alpha_k\rho_k e_k}{\alpha_k} - \Pi_k)/\Gamma_k \\
    \text{For new SGEOS, changed to } p_k = (\frac{\alpha_k\rho_k e_k - (\alpha_k\rho_k)q_k}{\alpha_k} - \Pi_k)/\Gamma_k\\
    p_k = \max(p_k, -(1-10^{-8})p^{\infty}_k + 10^{-8})
    $$

    else

    $$
    p_k = 0
    $$

    $p = \sum\alpha_k p_k$

    - interative process

    $$
    f(p) = 10^{-9} \\
    f'(p) = 10^9 \\
    \rho_k = 0
    $$

    loop when $f(p) > 10^{-10}$

    $$
    p = p - \frac{f(p)}{f'(p)} \\
    p = \max(p, -(1-10^{-8})p^{\infty}_k + 1) ?? \\
    f(p) = -1 \\
    f'(p) = 0 \\
    if \alpha_k > \epsilon, \rho_k = \frac{\alpha_k\rho_k}{\max(\alpha_k, \epsilon)}(\frac{p+p^{\infty}_k}{p_k+p^{\infty}_k})^{\frac{1}{\gamma_k}} \\
    f(p) = -1 + \sum\frac{\alpha_k\rho_k}{\rho_k} \\
    f'(p) = -\sum\frac{\alpha_k\rho_k}{\gamma_k\rho_k(p + p^{\infty}_k)}
    $$

    - update volume fraction

    $$
    if \alpha_k > \epsilon, \alpha_k = \frac{\alpha_k\rho_k}{\rho_k}
    $$

3. Mixture-total-energy correction

$$
\rho = \sum\alpha_k\rho_k \\
\Gamma = \sum\alpha_k\Gamma_k \\
\Pi = \sum\alpha_k\Pi_k \\
\text{For new SGEOS, add } mQ = \rho q = \sum(\alpha_k\rho_k)q_k \\
dyn_p = \frac{(\rho u)^2}{2\max(\rho,\epsilon)} \\
p = \frac{E - dyn_p - \Pi}{\Gamma} \\
\text{For new SGEOS, changed to } p = \frac{E - dyn_p - \rho q - \Pi}{\Gamma} \\
\alpha_k\rho_k e_k = \alpha_k(\Gamma_k p + \Pi_k) \\
\text{For new SGEOS, changed to } \alpha_k\rho_k e_k = \alpha_k(\Gamma_k p + \Pi_k) + (\alpha_k\rho_k)q_k
$$

## s_reconstruct_cell_boundary_values

to WENO-reconstruct the left and right cell-boundary values from the cell-averaged variables.

Variables:

- v_vf:
  - Cell-average variables
  - (sys_size) of scalar_field
  - IN
- vL_x, vR_x,vL_y,...:
  - left and right WENO-reconstructed cell-average variables
  - `(:,:,:,:)`
  - INOUT
- norm_dir:
  - splitting coordinate direction
- weno_dir:
  - WENO reconstruction direction

1. set direction

|     | norm_dir = 1   | 2 | 3 |
| ---:| :----: |:----:| :---   |
|  is1 | ix | iy | iz  |
|  is2 | iy | ix | iz  |
|  is3 | iz | iy | ix  |
|  weno_dir | 1 | 2 | 3 |
|  is1%beg | ix%beg + weno_polyn | iy%beg + weno_polyn | iz%beg + weno_polyn  |
|  is1%end | ix%end - weno_polyn | iy%end - weno_polyn | iz%end - weno_polyn  |

2. `s_weno`
   - v_vf
   - vL_x,vL_y,vL_z
   - vR_x,vR_y,vR_z

## s_finalize_rhs_module

deallocate variables and disassociate procedures.
