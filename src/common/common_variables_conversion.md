# m_variables_conversion.fpp

- conversion of the conservative variables into primitive ones and vice versa;
- obtain the mixture variables;
- compute pressure.

## s_compute_pressure (pressure calculation)

- model_eqns /= 4 and bubbles /= true

$$
dyn\_p = \frac{1}{2}\rho |u|^2 \\
energy = E = \rho e + \frac{1}{2}\rho |u|^2 \\
pi\_inf = \Pi = \frac{\gamma p^{\infty}}{\gamma - 1} \\
gamma = \Gamma = \frac{1}{\gamma - 1} \\
qv = \rho q
$$
$$
p = \frac{E - \frac{1}{2}\rho |u|^2 - \Pi - \rho q}{\Gamma}
$$

## s_convert_species_to_mixture_variables (mixture variables calculation)

$$
rho = \rho = \sum{\alpha_k\rho_k} \\
gamma = \Gamma = \sum{\alpha_k\Gamma_k} \\
pi\_inf = \Pi = \sum{\alpha_k\Pi_k} \\
qv = \rho q = \sum{\alpha_k\rho_k q_k}
$$

- mpp_lim
  - restrict $\alpha_k$ within [0,1] while keeping $\sum{\alpha_k}=1$;
  - keep $\alpha_k\rho_k >= 0$

## s_convert_conservative_to_primitive_variables

- input: qK_cons_vf
- output: qK_prim_vf

advxb: index of advecion fraction of phase 1.
contxe: index of the last continuity equation.

1. Use s_convert_species_to_mixture_variables_acc to compute mixture variables

2. Use s_compute_pressure to compute pressure

## s_convert_primitive_to_conservative_variables

- input: q_prim_vf
- output: q_cons_vf

1. Use s_convert_to_mixture_variables to copmute mixture variables

2. Compute energy from pressure

- model_eqns /= 4 and bubbles /= true

$$
E = \Gamma p + \Pi + \frac{1}{2}\rho |u|^2 + \rho q
$$

- model_eqns = 3

$$
\alpha_k\rho_k e_k = \alpha_k (\Gamma_k p + \Pi_k) + \alpha_k\rho_k q_k
$$

## s_convert_primitive_to_flux_variables

- input: qK_prim_vf
- output: FK_vf, FK_src_vf

1. Use s_convert_species_to_mixture_variables_acc to compute mixture variables

$$
E_K = \Gamma_K p_K +\Pi_K + \rho_K |u_K|^2/2 + (\rho q)_K
$$

2. Calculate flux

- mass flux

$$
F_K(\alpha_k\rho_k) = (\alpha_k\rho_k)_K u_K
$$

- momentom flux

$$
F_K(\rho u) = \rho_K u_K u_K + p_K
$$

- energy flux

$$
F_K(E) = u_K(E_K + p_K)
$$

- advection flux
  - riemann_solver == 1 (HLL)
    $$
    F_K(\alpha_k) = 0 \\
    F_{src,K}(\alpha_k) = (\alpha_k)_K
    $$
  - riemann_solver == 2 (HLLC)
    $$
    F_K(\alpha_k) = (\alpha_k)_K u_K \\
    F_{src,K}(\alpha_k) = u_K
    $$