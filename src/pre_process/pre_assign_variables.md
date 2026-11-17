# m_assign_variables.f90

Note: for the smoothing of the patch's boundaries, please read Shyue's work (1998)

## s_initialize_assign_variables_module

point s_assign_patch_primitive_variables to:

- Gamma/pi_inf model: s_assign_patch_mixture_primitive_variables
- Volume fraction model: s_assign_patch_species_primitive_variables

## s_assign_patch_mixture_primitive_variables (Gamma/pi_inf model)

...

## s_assign_patch_species_primitive_variables (Volume fraction model)

1. restore original primitive variables

    - orig_prim_vf
    - `s_convert_to_mixture_variables` to get
      - $\rho_{orig}$
      - $\Gamma_{orig}$
      - $\Pi_{orig}$
      - $q_{orig}$

2. compute current patch variables for smoothing

    - $\alpha_k\rho_k,\alpha_k$
    - `s_convert_to_mixture_variables` to get
      - $\rho_{patch}$
      - $\Gamma_{patch}$
      - $\Pi_{patch}$
      - $q_{patch}$

3. compute smoothing patch variables for smoothing

    - $\alpha_k\rho_k,\alpha_k$
    - `s_convert_to_mixture_variables` to get
      - $\rho_{smooth}$
      - $\Gamma_{smooth}$
      - $\Pi_{smooth}$
      - $q_{smooth}$

4. smoothing

    $$
    p = \eta p_{patch} + (1-\eta)p_{orig} \\
    \alpha_k = \eta (\alpha_k)_{patch} + (1-\eta)(\alpha_k)_{orig} \\
    \alpha_k\rho_k = \eta (\alpha_k\rho_k)_{patch} + (1-\eta)(\alpha_k\rho_k)_{orig} \\
    u = \eta u_{patch} + (1-\eta) v_{orig}
    $$

    - `s_convert_to_mixture_variables` to get
      - $\rho$
      - $\Gamma$
      - $\Pi$
      - $q$

    - for six-equation model (model_eqns = 3)
      - ($\alpha_1\rho_1 e_1$): $p_1 = p$ 
      - ($\alpha_2\rho_2 e_2$): $p_2 = p$

## s_finalize_assign_varaibles_module

point `s_assign_patch_primitive_variables` => null()