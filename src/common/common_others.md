# Other simple functions

## m_constants.fpp (constant parameters for utilization)

dflt_real = $-10^6$ <br>
sgm_eps = $10^{-16}$ <br>
small_alf = $10^{-11}$ <br>
verysmall = $10^{-12}$

## m_derived_types.fpp (custom-defined types)

- field_position: $(:,:,:)$
- scalar_field: $(:,:,:)$
- integer_field: $(:,:,:)$, integer
- vector_field: $(:)$, scalar_field
- physical_parameters:
  - gamma   $\frac{1}{\gamma - 1}$
  - pi_inf  $\frac{\gamma p^{\infty}}{\gamma - 1}$
  - Re(2)   Reynolds number
  - cv      $C_v$
  - qv      $q$
  - qvp     $q^{\prime}$
  - mul0    $\mu_0$
  - ss      bubble surface tension
  - pv      bubble vapour pressure
  - gamma_v bubble constants
  - M_v     bubble constants
  - mu_v    bubble constants
  - k_v     bubble constants
  - G

## m_helper.fpp (utilized in subgrid bubble model)

### s_compute_finite_difference_coefficients

Variables

- s_cc: locations of the cell-centers in the s-coordinate direction
- fd_coeff_s: finite-differece coefficients in the s-coordinate direction

1. fd_order_in = 1:

| Term     | Coefficient              |
|----------|--------------------------|
| `0`| `-1/Δx`               |
| `1`| `1/Δx`              |

2. order = 2

| Term     | Coefficient              |
|----------|--------------------------|
| `-1`| `-1/(2Δx)`              |
| `0` | `0`                      |
| `1` | `1/(2Δx)`               |

3. order = 4

| Term     | Coefficient              |
|----------|--------------------------|
| `-2`| `1/(12Δx)`               |
| `-1`| `-8/(12Δx)`              |
| `0` | `0`                      |
| `1` | `8/(12Δx)`               |
| `2` | `-1/(12Δx)`              |

### s_comp_n_from_prim (bubble number density computation from primitive variables)

...