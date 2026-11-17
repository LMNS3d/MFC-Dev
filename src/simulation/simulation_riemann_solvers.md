# m_riemann_solvers.fpp

Variables:

- flux_rsx_vf, flux_rsy_vf, flux_rsz_vf
  - `(:,:,:,:)`
  - dimension-splitting cell-boundary fluxes
- flux_src_rsx_vf, flux_src_rsy_vf, flux_src_rsz_vf
  - `(:,:,:,:)`
  - dimension-splitting cell-boudary sources
- vel_src_rsx_vf, vel_src_rsy_vf, vel_src_rsz_vf
  - `(:,:,:,:)`
  - dimension-splitting cell-boudary velocity

## s_hll_riemann_solver

Variables:

- qL_prim_rsx_vf, qR_prim_rsx_vf, qL_prim_rsy_vf, ...
  - `(startx:,starty:,startz:,1:)`, first three for spatial indexes, last for equation
  - INOUT
- q_prim_vf
  - (sys_size) of scalar_field
  - IN
- qL_prim_vf, qR_prim_vf
  - `(:)` of scalar_field
  - INOUT
- flux_vf, flux_src_vf
  - (sys_size) of scalar_field
  - INOUT

1. `s_populate_riemann_states_variables_buffers`
    - populate buffers of the left and right Riemann states variables based on boundary conditions

2. `s_initialize_riemann_solver`
    - reshaping inputted data based on dimension-splitting direction

3. Riemann solver by dimensional splitting for cell-boundary fluxes and sources

    Note:

    - qL_prim_rsx_vf(j, k, l, i) corresponds to $prim_{j+\frac{1}{2}}^L$
    - qR_prim_rsx_vf(j+1, k, l, i) corresponds to $prim_{j+\frac{1}{2}}^R$

    1. load reconstructed variables
     $(\alpha_k\rho_k)_L, (\alpha_k\rho_k)_R, V_L, V_R, (\alpha_k)_L, (\alpha_k)_R, p_L, p_R$

      if mpp_lim:
      $$
      (\alpha_k\rho_k)_L = \max(0,(\alpha_k\rho_k)_L) \\
      (\alpha_k\rho_k)_R = \max(0,(\alpha_k\rho_k)_R) \\
      (\alpha_k)_L = \min(\max(0,(\alpha_k)_L),1) \\
      (\alpha_k)_R = \min(\max(0,(\alpha_k)_R),1) \\
      (\alpha_k)_L = \frac{(\alpha_k)_L}{\max(\epsilon,\sum(\alpha_k)_L)} \\
      (\alpha_k)_R = \frac{(\alpha_k)_R}{\max(\epsilon,\sum(\alpha_k)_R)} \\
      $$

    2. calculate mixing left and right cell-boundary values

        $$
        \rho_L = \sum(\alpha_k\rho_k)_L, \rho_R = \sum(\alpha_k\rho_k)_R \\
        \Gamma_L = \sum(\alpha_k)_L\Gamma_k, \Gamma_R = \sum(\alpha_k)_R\Gamma_k \\
        \Pi_L = \sum(\alpha_k)_L\Pi_k, \Pi_R = \sum(\alpha_k)_R\Pi_k \\
        Q_L = \sum(\alpha_k\rho_k)_Lq_k, Q_R = \sum(\alpha_k\rho_k)_Rq_k \\
        (V_{rms}^2)_L = |V_L|^2, (V_{rms}^2)_R = |V_R|^2
        $$

        $$
        E_L = \Gamma_Lp_L + \Pi_L + \frac{1}{2}\rho_L|V_L|^2 + Q_L \\
        E_R = \Gamma_Rp_R + \Pi_R + \frac{1}{2}\rho_R|V_R|^2 + Q_R \\
        H_L = (E_L + p_L)/\rho_L \\
        H_R = (E_R + p_R)/\rho_R
        $$

    3. calculate wave speeds

     - `compute_average_state`
       - to calculate average values of the left and right-biased variables
     - `s_compute_speed_of_sound`
       - to get $c_L,c_R,c_{avg}$

        if wave_speeds = 1 (direct wave speed)

        $$
        s_L = \min(u_L - c_L, u_R - c_R) \\
        s_R = \max(u_R + c_R, u_L + c_L) \\
        s^{\star} = \frac{p_R - p_L + \rho_Lu_L(s_L - u_L) - \rho_Ru_R(s_R - u_R)}{\rho_L(s_L - u_L) - \rho_R(s_R - u_R)}
        $$

        if wave_speeds = 2

        $$
        p^{\star}_L = \frac{p_L + p_R + \rho_{avg}c_{avg}(u_L - u_R)}{2} \\
        p^{\star}_R = p^{\star}_L \\
        Ms_L = \max\left(1, \sqrt{1 + \frac{\frac{1}{2} + \Gamma_L}{1 + \Gamma_L}\left(\frac{p^{\star}_L - p_L}{p_L + \frac{\Pi_L}{1 + \Gamma_L}}\right)}\right) \\
        Ms_R = \max\left(1, \sqrt{1 + \frac{\frac{1}{2} + \Gamma_R}{1 + \Gamma_R}\left(\frac{p^{\star}_R - p_R}{p_R + \frac{\Pi_R}{1 + \Gamma_R}}\right)}\right) \\
        s_L = u_L - c_L\cdot Ms_L \\
        s_R = u_R + c_R\cdot Ms_R \\
        s^{\star} = \frac{u_L + u_R}{2} + \frac{p_L - p_ R}{2 \rho_{avg} c_{avg}}
        $$

        $$
        s_- = \min(0,s_L), s_+ = \max(0,s_R) \\
        \xi_- = \left(\frac{1}{2} + \text{sign}(\frac{1}{2}, s_L)\right) + \\
        \left(\frac{1}{2} - \text{sign}(\frac{1}{2}, s_L)\right) \times \left(\frac{1}{2} + \text{sign}(\frac{1}{2}, s_R)\right) \\
        if s_R>0,\xi_- = 1 \\
        \xi_+ = \left(\frac{1}{2} - \text{sign}(\frac{1}{2}, s_R)\right) + \\
        \left(\frac{1}{2} - \text{sign}(\frac{1}{2}, s_L)\right) \times \left(\frac{1}{2} + \text{sign}(\frac{1}{2}, s_R)\right) \\
        if s_L<0,\xi_+ = 1
        $$

    4. calculate fluxes and sources

$$
F_{i+\frac{1}{2}}(\alpha_k\rho_k) = \frac{s_-(\alpha_k\rho_k)_Ru_R - s_+(\alpha_k\rho_k)_Lu_L + s_-s_+((\alpha_k\rho_k)_L - (\alpha_k\rho_k)_R)}{s_- - s_+} \\
F_{i+\frac{1}{2}}(\rho V_j) = \frac{s_-(\rho_RV_{sp,R}V_{j,R} + p_R\delta_{sp,j}) - s_+(\rho_LV_{sp,L}V_{j,L} + p_L\delta_{sp,j})+s_-s_+(\rho_LV_{j,L} - \rho_RV_{j,R})}{s_- - s_+} \\
F_{i+\frac{1}{2}}(E) = \frac{s_-u_R(E_R + p_R) - s_+u_L(E_L + p_L) + s_-s_+(E_L -E_R)}{s_- - s_+} \\
F_{i+\frac{1}{2}}(\alpha_k) = \frac{s_-s_+((\alpha_k)_L - (\alpha_k)_R)}{s_- - s_+} \\
S_{i+\frac{1}{2}}(\alpha_k) = \frac{s_-(\alpha_k)_R - s_+(\alpha_k)_L}{s_- - s_+} \\
V_{j,i+\frac{1}{2}} = \frac{\xi_-(\rho_LV_{j,L}(s_L - V_{sp,L}) - p_L\delta_{sp,j}) - \xi_+(\rho_RV_{j,R}(s_R - V_{sp,R}) - p_R\delta_{sp,j})}{\xi_-\rho_L(s_L - V_{sp,L}) - \xi_+\rho_R(s_R - V_{sp,R})}
$$

4. `s_finalize_riemann_solver`

Reshape flux and source variables.

- norm_dir = 2
  - flux_vf(i)%sf(k, j, l) = flux_rsy_vf(j, k, l, i)
  - flux_src_vf(advxb:advxe)%sf(k, j, l) = flux_src_rsy_vf(j, k, l, advxb:advxe)
- norm_dir = 3
  - flux_vf(i)%sf(l, k, j) = flux_rsy_vf(j, k, l, i)
  - flux_src_vf(advxb:advxe)%sf(l, k, j) = flux_src_rsy_vf(j, k, l, advxb:advxe)
- norm_dir = 1
  - flux_vf(i)%sf(j, k, l) = flux_rsy_vf(j, k, l, i)
  - flux_src_vf(advxb:advxe)%sf(j, k, l) = flux_src_rsy_vf(j, k, l, advxb:advxe)

## s_hllc_riemann_solver

Note: other steps are similar to HLL, here we focus on dimensional splitting parts.

- model_eqns = 3 (six-equation model)

  1. load reconstructed variables
       $(\alpha_k\rho_k)_L, (\alpha_k\rho_k)_R, V_L, V_R, (\alpha_k)_L, (\alpha_k)_R, p_L, p_R$

        $$
        V_{rms,L}^2 = \sum V_L^2 \\
        V_{rms,R}^2 = \sum V_R^2
        $$

  2. calculate mixing left and right cell-boundary values
  3. calculate wave speeds

     - `compute_average_state`
       - to calculate average values of the left and right-biased variables
     - `s_compute_speed_of_sound`
       - to get $c_L,c_R,c_{avg}$

        if wave_speeds = 1 (direct wave speed)

        $$
        s_L = \min(u_L - c_L, u_R - c_R) \\
        s_R = \max(u_R + c_R, u_L + c_L) \\
        s^{\star} = \frac{p_R - p_L + \rho_Lu_L(s_L - u_L) - \rho_Ru_R(s_R - u_R)}{\rho_L(s_L - u_L) - \rho_R(s_R - u_R)}
        $$

        if wave_speeds = 2, ...

  4. calculate fluxes based on wave speeds

        Note: <br>
        norm_dir = 1: $V_{sp} = V_1 = u$; <br>
        norm_dir = 2: $V_{sp} = V_2 = v$; <br>
        norm_dir = 3: $V_{sp} = V_3 = w$.

        - $s_L \ge 0$
        $$
        p^{\star} = p_L \\
        F_{i+\frac{1}{2}}(\alpha_k) = (\alpha_k)_Ls^{\star} \\
        F_{i+\frac{1}{2}}(\alpha_k\rho_k) = (\alpha_k\rho_k)_Lu_L \\
        F_{i+\frac{1}{2}}(\alpha_k\rho_k e_k) = ((\alpha_k)_L(\Gamma_kp_L + \Pi_k) + (\alpha_k\rho_k)_Lq_k)u_L \\
        F_{i+\frac{1}{2}}(\rho V_j) = \rho_LV_{sp,L}V_{j,L}+p_L\delta_{sp,j} \\
        V_{j,i+\frac{1}{2}} = V_{j,L} + \delta_{sp,j}(s^{\star}-V_{j,L}) \\
        F_{i+\frac{1}{2}}(E) = (E_L + p_L)u_L
        $$
        - $s_R \le 0$
        $$
        p^{\star} = p_R \\
        F_{i+\frac{1}{2}}(\alpha_k) = (\alpha_k)_Rs^{\star} \\
        F_{i+\frac{1}{2}}(\alpha_k\rho_k) = (\alpha_k\rho_k)_Ru_R \\
        F_{i+\frac{1}{2}}(\alpha_k\rho_k e_k) = ((\alpha_k)_R(\Gamma_kp_R + \Pi_k) + (\alpha_k\rho_k)_Rq_k)u_R \\
        F_{i+\frac{1}{2}}(\rho V_j) = \rho_RV_{sp,R}V_{j,R}+p_R\delta_{sp,j} \\
        V_{j,i+\frac{1}{2}} = V_{j,R} + \delta_{sp,j}(s^{\star}-V_{j,R}) \\
        F_{i+\frac{1}{2}}(E) = (E_R + p_R)u_R
        $$
        - $s_L < 0 \le s^{\star}$
        $$
        \xi_L = \frac{s_L - u_L}{s_L - s^{\star}} \\
        \rho^{\star} = \rho_L\xi_L \\
        E^{\star} = \xi_L(E_L + (s^{\star} - u_L)(\rho_Ls^{\star} + \frac{p_L}{s_L - u_L})) \\
        p^{\star} = \rho_L(s_L - u_L)(s^{\star}-u_L) + p_L \\
        p^{\star}_k = (p_L + \frac{\Pi_k}{1+\Gamma_k})\xi_L^{\frac{1}{\Gamma_k}+1} - \frac{\Pi_k}{1+\Gamma_k} \\
        F_{i+\frac{1}{2}}(\alpha_k) = (\alpha_k)_Ls^{\star} \\
        F_{i+\frac{1}{2}}(\alpha_k\rho_k) = (\alpha_k\rho_k)_L\xi_Ls^{\star} \\
        F_{i+\frac{1}{2}}(\alpha_k\rho_k e_k) = ((\alpha_k)_L(\Gamma_kp_k^{\star} + \Pi_k) + (\alpha_k\rho_k)_Lq_k)s^{\star} \\
        F_{i+\frac{1}{2}}(\rho V_j) = \rho^{\star}s^{\star}(s^{\star}\delta_{sp,j} + V_{j,L}(1-\delta_{sp,j}))+p^{\star}\delta_{sp,j} \\
        V_{j,i+\frac{1}{2}} = V_{j,L} + \delta_{sp,j}(s^{\star}\xi_L-V_{j,L}) \\
        F_{i+\frac{1}{2}}(E) = (E^{\star} + p^{\star})u^{\star}
        $$
        - $s^{\star} < 0 < s_R$
        $$
        \xi_R = \frac{s_R - u_R}{s_R - s^{\star}} \\
        \rho^{\star} = \rho_R \xi_R \\
        E^{\star} = \xi_R(E_R + (s^{\star} - u_R)(\rho_Rs^{\star} + \frac{p_R}{s_R - u_R})) \\
        p^{\star} = \rho_R(s_R - u_R)(s^{\star}-u_R) + p_R \\
        p^{\star}_k = (p_R + \frac{\Pi_k}{1+\Gamma_k})\xi_R^{\frac{1}{\Gamma_k}+1} - \frac{\Pi_k}{1+\Gamma_k} \\
        F_{i+\frac{1}{2}}(\alpha_k) = (\alpha_k)_Rs^{\star} \\
        F_{i+\frac{1}{2}}(\alpha_k\rho_k) = (\alpha_k\rho_k)_R\xi_Rs^{\star} \\
        F_{i+\frac{1}{2}}(\alpha_k\rho_k e_k) = ((\alpha_k)_R(\Gamma_kp_k^{\star} + \Pi_k) + (\alpha_k\rho_k)_Rq_k)s^{\star} \\
        F_{i+\frac{1}{2}}(\rho V_j) = \rho^{\star}s^{\star}(s^{\star}\delta_{sp,j} + V_{j,R}(1-\delta_{sp,j})) + p^{\star}\delta_{sp,j} \\
        V_{j,i+\frac{1}{2}} = V_{j,R} + \delta_{sp,j}(s^{\star}\xi_R-V_{j,R}) \\
        F_{i+\frac{1}{2}}(E) = (E^{\star} + p^{\star})u^{\star}
        $$

  5. calculate source term

        flux_src_rsx/y/z = $V_{sp,i+\frac{1}{2}}$

- model_eqns = 2 (five-equation model)

  1. load reconstructed variables
       $(\alpha_k\rho_k)_L, (\alpha_k\rho_k)_R, V_L, V_R, (\alpha_k)_L, (\alpha_k)_R, p_L, p_R$
  2. calculate mixing left and right cell-boundary values
  3. calculate wave speeds

     - `compute_average_state`
       - to calculate average values of the left and right-biased variables
     - `s_compute_speed_of_sound`
       - to get $c_L,c_R,c_{avg}$

        if wave_speeds = 1 (direct wave speed)

        $$
        s_L = \min(u_L - c_L, u_R - c_R) \\
        s_R = \max(u_R + c_R, u_L + c_L) \\
        s^{\star} = \frac{p_R - p_L + \rho_Lu_L(s_L - u_L) - \rho_Ru_R(s_R - u_R)}{\rho_L(s_L - u_L) - \rho_R(s_R - u_R)}
        $$

        if wave_speeds = 2, ...

    $$
    s_- = \min(0,s_L). s_+ = \max(0,s_R) \\
    \xi_L = \frac{s_L - u_L}{s_L - s^{\star}}, \xi_R = \frac{s_R - u_R}{s_R - s^{\star}} \\
    \xi_- = \frac{1}{2} + sign(\frac{1}{2},s^{\star}) \\
    \xi_+ = \frac{1}{2} - sign(\frac{1}{2},s^{\star})
    $$

  4. calculate fluxes based on wave speeds

        $$
        \begin{align*}
            F_{i+\frac{1}{2}}(\alpha_k\rho_k) &= \xi_-(\alpha_k\rho_k)_L(u_L + s_-(\xi_L-1))  \\
            &+ \xi_+(\alpha_k\rho_k)_R(u_R + s_+(\xi_R - 1))
        \end{align*} \\
        \begin{align*}
            F_{i+\frac{1}{2}}(\rho V_j) &= \xi_-(\rho_L(V_{sp,L}V_{j,L} \\
            & + s_-(\xi_L(s^{\star}\delta_{sp,j} + V_{j,L}(1-\delta_{sp,j})) - V_{j,L})) + p_L\delta_{sp,j}) \\
            & + \xi_+(\rho_R(V_{sp,R}V_{j,R} \\
            & + s_+(\xi_R(s^{\star}\delta_{sp,j} + V_{j,R}(1-\delta_{sp,j})) - V_{j,R})) + p_R\delta_{sp,j})
        \end{align*} \\
        \begin{align*}
            F_{i+\frac{1}{2}}(E) & = \xi_-(u_L(E_L + p_L) \\
            & + s_-(\xi_L(E_L + (s^{\star}-u_L)(\rho_Ls^{\star} + \frac{p_L}{s_L - u_L})) - E_L)) \\
            & + \xi_+(u_R(E_R + p_R) \\
            & + s_+(\xi_R(E_R + (s^{\star} - u_L)(\rho_R s^{\star} + \frac{p_R}{s_R - u_R})) - E_R))
        \end{align*} \\
        \begin{align*}
            F_{i+\frac{1}{2}}(\alpha_k) & = \xi_-(\alpha_k)_L(u_L + s_-(\xi_L - 1)) \\
            & + \xi_+(\alpha_k)_R(u_R + s_+(\xi_R - 1))
        \end{align*} \\
        V_{j,i+\frac{1}{2}} = \xi_-(V_{j,L} +\delta_{sp,j}s_-(\xi_L - 1)) + \xi_+(V_{j,R} + \delta_{sp,j}s_+(\xi_R - 1))
        $$

  5. calculate source

        flux_src_rs = $V_{sp,i+\frac{1}{2}}$

## s_initialize_riemann_solvers_module

allocation of memory, association of pointers, computation of parameters.

point `s_riemann_solver` to:

- riemann_solver = 1: `s_hll_riemann_solver`
- riemann_solver = 2: `s_hllc_riemann_solver`

Allocation:

- `flux_rsx/y/z_vf(is1,is2,is3,1:sys_size)`
- `flux_src_rsx/y/z_vf(is1,is2,is3,advxb:sys_size)`
- `vel_src_rsx/y/z_vf(is1,is2,is3,1:num_dims)`

for rsx

|   | is1  | is2  | is3  |
|------|------|------|------|
| beg | -1 | 0 | 0 |
| end | m | n | p |

for rsy

|   | is1  | is2  | is3  |
|------|------|------|------|
| beg | -1 | 0 | 0 |
| end | n | m | p |

for rsz

|   | is1  | is2  | is3  |
|------|------|------|------|
| beg | -1 | 0 | 0 |
| end | p | n | m |

## s_populate_riemann_states_variables_buffers

populate the buffers of the left and right Riemann states variables, depending on the boundary conditions.

|  (dir_idx, dir_flg, is) | 1  | 2  | 3  |
|------|------|------|------|
| norm_dir=1 | `(1,1,ix)` | `(2,0,iy)` | `(3,0,iz)` |
| norm_dir=2 | `(2,0,iy)` | `(1,1,ix)` | `(3,0,iz)` |
| norm_dir=3 | `(3,0,iz)` | `(1,0,iy)` | `(2,1,ix)` |

## s_initialize_riemann_solver

Nothing utilized.

## s_finalize_riemann_solver

- norm_dir = 2
  - flux_vf(i)%sf(k, j, l) = flux_rsy_vf(j, k, l, i)
  - flux_src_vf(`advxb`)%sf(k, j, l) = flux_src_rsy_vf(j, k, l, advxb)
  - if `HLL`, flux_src_vf(`advxb:advxe`)%sf(k, j, l) = flux_src_rsy_vf(j, k, l, advxb:advxe)
- norm_dir = 3
  - flux_vf(i)%sf(l, k, j) = flux_rsy_vf(j, k, l, i)
  - flux_src_vf(`advxb`)%sf(l, k, j) = flux_src_rsy_vf(j, k, l, advxb)
  - if `HLL`, flux_src_vf(`advxb:advxe`)%sf(l, k, j) = flux_src_rsy_vf(j, k, l, advxb:advxe)
- norm_dir = 1
  - flux_vf(i)%sf(j, k, l) = flux_rsy_vf(j, k, l, i)
  - flux_src_vf(`advxb`)%sf(j, k, l) = flux_src_rsy_vf(j, k, l, advxb)
  - if `HLL`, flux_src_vf(`advxb:advxe`)%sf(j, k, l) = flux_src_rsy_vf(j, k, l, advxb:advxe)

## s_finalize_riemann_solvers_module

deallocation and disassociation.

`s_riemann_solver` => null()

deallocate vel_src_rsx/y/z_vf, flux_rsx/y/z_vf, flux_src_rsx/y/z_vf
