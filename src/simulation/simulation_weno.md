# m_weno.fpp

Variables

- v_rs_ws_x, v_rs_ws_y, v_rs_ws_z
  - `(:,:,:,:)`
  - reshaped cell-average variables for WENO reconstruction
- poly_coef_cbL_x, poly_coef_cbL_y, poly_coef_cbL_z
  - `(:,:,:)`: second identifies the polynomial, third identifies the position of its coefficients, first denotes the cell-location
  - polynomial coefficients at the left cell-boundaries
- poly_coef_cbR_x, poly_coef_cbR_y, poly_coef_cbR_z
  - `(:,:,:)`: second identifies the polynomial, third identifies the position of its coefficients, first denotes the cell-location
  - polynomial coefficients at the right cell-boundaries
- d_cbL_x, d_cbL_y, d_cbL_z
  - `(:,:)`: first identifies the weight, last denotes the cell-location
  - ideal weights at the left cell-boundaries
- d_cbR_x, d_cbR_y, d_cbR_z
  - `(:,:)`: first identifies the weight, last denotes the cell-location
  - ideal weights at the right cell-boundaries
- beta_coef_x, beta_coef_y, beta_coef_z
  - `(:,:,:)`: first identifies the smoothness indicator, second identifies the position of coefficients, last denotes the cell-location
  - smoothness indicator coefficients

## s_initialize_weno_module

memory allocation, pointers association, parameters computation.

|   | is1_weno  | is2_weno  | is3_weno  |
|------|------|------|------|
| beg | -buff_size | -buff_size | -buff_size |
| end | m + buff_size | n + buff_size | p + buff_size |

allocate

|   | poly_coef_cbL_x  | d_cbL_x  | beta_coef_x  |v_rs_ws_x |
|------|------|------|------| ------ |
| 1 | is1_weno $\pm$ weno_polyn | 0:weno_polyn | is1_weno $\pm$ weno_polyn | is1_weno |
| 2 | 0:weno_polyn | is1_weno $\pm$ weno_polyn | 0:weno_polyn | is2_weno |
| 3 | 0:weno_polyn - 1 |  | 0:2*(weno_polyn - 1) | is3_weno |
| 4 |  |  |  | 1:sys_size |

call `s_compute_weno_coefficients`

## s_compute_weno_coefficients

compute the grid dependent coefficients of the WENO reconstruction

- weno_order = 3

$$
\text{poly\_coef\_cbR\_}{XYZ}(i + 1, 0, 0) = \frac{s_{\text{cb}}(i) - s_{\text{cb}}(i + 1)}{s_{\text{cb}}(i) - s_{\text{cb}}(i + 2)} = \frac{1}{2} \\
\text{poly\_coef\_cbR\_}{XYZ}(i + 1, 1, 0) = \frac{s_{\text{cb}}(i) - s_{\text{cb}}(i + 1)}{s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 1)} = \frac{1}{2}
$$

$$
\text{poly\_coef\_cbL\_}{XYZ}(i + 1, 0, 0) = -\frac{s_{\text{cb}}(i) - s_{\text{cb}}(i + 1)}{s_{\text{cb}}(i) - s_{\text{cb}}(i + 2)} = -\frac{1}{2} \\
\text{poly\_coef\_cbL\_}{XYZ}(i + 1, 1, 0) = -\frac{s_{\text{cb}}(i) - s_{\text{cb}}(i + 1)}{s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 1)} = -\frac{1}{2}
$$

$$
d_{\text{cbR\_}{XYZ}}(0, i + 1) = \frac{s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 1)}{s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 2)} = \frac{2}{3} \\
d_{\text{cbL\_}{XYZ}}(0, i + 1) = \frac{s_{\text{cb}}(i - 1) - s_{\text{cb}}(i)}{s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 2)} = \frac{1}{3} \\
d_{\text{cbR\_}{XYZ}}(1, i + 1) = 1 - d_{\text{cbR\_}{XYZ}}(0, i + 1) = \frac{1}{3} \\
d_{\text{cbL\_}{XYZ}}(1, i + 1) = 1 - d_{\text{cbL\_}{XYZ}}(0, i + 1) = \frac{2}{3}
$$

$$
\beta_{\text{coef\_}{XYZ}}(i + 1, 0, 0) = 4 \left(\frac{s_{\text{cb}}(i) - s_{\text{cb}}(i + 1)}{s_{\text{cb}}(i) - s_{\text{cb}}(i + 2)}\right)^2 = 1 \\
\beta_{\text{coef\_}{XYZ}}(i + 1, 1, 0) = 4 \left(\frac{s_{\text{cb}}(i) - s_{\text{cb}}(i + 1)}{s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 1)}\right)^2 = 1
$$

- weno_order = 5

$$
\text{poly\_coef\_cbR\_${XYZ}$} (i + 1, 0, 0)  \\
= \frac{(s_{\text{cb}}(i) - s_{\text{cb}}(i + 1)) \cdot (s_{\text{cb}}(i + 1) - s_{\text{cb}}(i + 2))}{(s_{\text{cb}}(i) - s_{\text{cb}}(i + 3)) \cdot (s_{\text{cb}}(i + 3) - s_{\text{cb}}(i + 1))} = -\frac{1}{6}
$$

$$
\text{poly\_coef\_cbR\_${XYZ}$} (i + 1, 1, 0) \\
= \frac{(s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 1)) \cdot (s_{\text{cb}}(i + 1) - s_{\text{cb}}(i))}{(s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 2)) \cdot (s_{\text{cb}}(i + 2) - s_{\text{cb}}(i))} = \frac{2}{6}
$$

$$
\text{poly\_coef\_cbR\_${XYZ}$} (i + 1, 1, 1) \\
= \frac{(s_{\text{cb}}(i) - s_{\text{cb}}(i + 1)) \cdot (s_{\text{cb}}(i + 1) - s_{\text{cb}}(i + 2))}{(s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 1)) \cdot (s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 2))} = \frac{1}{6}
$$

$$
\text{poly\_coef\_cbR\_${XYZ}$} (i + 1, 2, 1) \\
= \frac{(s_{\text{cb}}(i) - s_{\text{cb}}(i + 1)) \cdot (s_{\text{cb}}(i + 1) - s_{\text{cb}}(i - 1))}{(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i)) \cdot (s_{\text{cb}}(i - 2) - s_{\text{cb}}(i + 1))} = -\frac{2}{6}
$$

$$
\text{poly\_coef\_cbR\_${XYZ}$} (i + 1, 0, 1) \\
= \frac{(s_{\text{cb}}(i) - s_{\text{cb}}(i + 2)) + (s_{\text{cb}}(i + 1) - s_{\text{cb}}(i + 3))}{(s_{\text{cb}}(i) - s_{\text{cb}}(i + 2)) \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i + 3))} \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i + 1)) = \frac{4}{6}
$$

$$
\text{poly\_coef\_cbR\_${XYZ}$} (i + 1, 2, 0) \\
= \frac{(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i + 1)) + (s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 1))}{(s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 1)) \cdot (s_{\text{cb}}(i + 1) - s_{\text{cb}}(i - 2))} \cdot (s_{\text{cb}}(i + 1) - s_{\text{cb}}(i)) = \frac{5}{6}
$$

$$
\text{poly\_coef\_cbL\_${XYZ}$} (i + 1, 0, 0) \\
= \frac{(s_{\text{cb}}(i + 1) - s_{\text{cb}}(i))(s_{\text{cb}}(i) - s_{\text{cb}}(i + 2))}{(s_{\text{cb}}(i) - s_{\text{cb}}(i + 3))(s_{\text{cb}}(i + 3) - s_{\text{cb}}(i + 1))} = \frac{2}{6}
$$

$$
\text{poly\_coef\_cbL\_${XYZ}$} (i + 1, 1, 0) \\
= \frac{(s_{\text{cb}}(i) - s_{\text{cb}}(i - 1))(s_{\text{cb}}(i) - s_{\text{cb}}(i + 1))}{(s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 2))(s_{\text{cb}}(i) - s_{\text{cb}}(i + 2))} = -\frac{1}{6}
$$

$$
\text{poly\_coef\_cbL\_${XYZ}$} (i + 1, 1, 1) \\
= \frac{(s_{\text{cb}}(i + 1) - s_{\text{cb}}(i))(s_{\text{cb}}(i) - s_{\text{cb}}(i + 2))}{(s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 1))(s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 2))} = -\frac{2}{6}
$$

$$
\text{poly\_coef\_cbL\_${XYZ}$} (i + 1, 2, 1) \\
= \frac{(s_{\text{cb}}(i - 1) - s_{\text{cb}}(i))(s_{\text{cb}}(i) - s_{\text{cb}}(i + 1))}{(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i))(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i + 1))} = \frac{1}{6}
$$

$$
\text{poly\_coef\_cbL\_${XYZ}$} (i + 1, 0, 1) \\
= \frac{(s_{\text{cb}}(i) - s_{\text{cb}}(i + 2)) + (s_{\text{cb}}(i) - s_{\text{cb}}(i + 3))}{(s_{\text{cb}}(i) - s_{\text{cb}}(i + 2))(s_{\text{cb}}(i) - s_{\text{cb}}(i + 3))} \cdot (s_{\text{cb}}(i + 1) - s_{\text{cb}}(i)) = -\frac{5}{6}
$$

$$
\text{poly\_coef\_cbL\_${XYZ}$} (i + 1, 2, 0) \\
= \frac{(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i)) + (s_{\text{cb}}(i - 1) - s_{\text{cb}}(i + 1))}{(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i + 1))(s_{\text{cb}}(i + 1) - s_{\text{cb}}(i - 1))} \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i + 1)) = -\frac{4}{6}
$$

$$
\text{d\_cbR\_${XYZ}$} (0, i + 1) = \frac{(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i + 1))(s_{\text{cb}}(i + 1) - s_{\text{cb}}(i - 1))}{(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i + 3))(s_{\text{cb}}(i + 3) - s_{\text{cb}}(i - 1))} = \frac{6}{20}
$$

$$
\text{d\_cbR\_${XYZ}$} (2, i + 1) = \frac{(s_{\text{cb}}(i + 1) - s_{\text{cb}}(i + 2))(s_{\text{cb}}(i + 1) - s_{\text{cb}}(i + 3))}{(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i + 2))(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i + 3))} = \frac{2}{20}
$$

$$
\text{d\_cbL\_${XYZ}$} (0, i + 1) = \frac{(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i))(s_{\text{cb}}(i) - s_{\text{cb}}(i - 1))}{(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i + 3))(s_{\text{cb}}(i + 3) - s_{\text{cb}}(i - 1))} = \frac{2}{20}
$$

$$
\text{d\_cbL\_${XYZ}$} (2, i + 1) = \frac{(s_{\text{cb}}(i) - s_{\text{cb}}(i + 2))(s_{\text{cb}}(i) - s_{\text{cb}}(i + 3))}{(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i + 2))(s_{\text{cb}}(i - 2) - s_{\text{cb}}(i + 3))} = \frac{6}{20}
$$

$$
\text{d\_cbR\_${XYZ}$} (1, i + 1) = 1.0 - \text{d\_cbR\_${XYZ}$} (0, i + 1) - \text{d\_cbR\_${XYZ}$} (2, i + 1) = \frac{12}{20}
$$

$$
\text{d\_cbL\_${XYZ}$} (1, i + 1) = 1.0 - \text{d\_cbL\_${XYZ}$} (0, i + 1) - \text{d\_cbL\_${XYZ}$} (2, i + 1) = \frac{12}{20}
$$

$$
\beta_{\text{coef\_}\text{${XYZ}$}}(i+1, 0, 0) \\
= \frac{4 \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i+1))^2 \cdot \left(10 \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i))^2 + (s_{\text{cb}}(i+1) - s_{\text{cb}}(i)) \cdot (s_{\text{cb}}(i+2) - s_{\text{cb}}(i+1)) + (s_{\text{cb}}(i+2) - s_{\text{cb}}(i+1))^2\right)}{(s_{\text{cb}}(i) - s_{\text{cb}}(i+3))^2 \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i+3))^2} = \frac{48}{36}
$$

$$
\beta_{\text{coef\_}\text{${XYZ}$}}(i+1, 0, 1) \\
= \frac{4 \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i+1))^2 \cdot \left(19 \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i))^2 - (s_{\text{cb}}(i+1) - s_{\text{cb}}(i)) \cdot (s_{\text{cb}}(i+3) - s_{\text{cb}}(i+1)) + 2 \cdot (s_{\text{cb}}(i+2) - s_{\text{cb}}(i)) \cdot ((s_{\text{cb}}(i+2) - s_{\text{cb}}(i)) + (s_{\text{cb}}(i+3) - s_{\text{cb}}(i+1)))\right)}{(s_{\text{cb}}(i) - s_{\text{cb}}(i+2)) \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i+3))^2 \cdot (s_{\text{cb}}(i+3) - s_{\text{cb}}(i+1))} = -\frac{132}{36}
$$

$$
\beta_{\text{coef\_}\text{${XYZ}$}}(i+1, 0, 2) \\
= \frac{4 \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i+1))^2 \cdot \left(10 \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i))^2 + (s_{\text{cb}}(i+1) - s_{\text{cb}}(i)) \cdot ((s_{\text{cb}}(i+2) - s_{\text{cb}}(i)) + (s_{\text{cb}}(i+3) - s_{\text{cb}}(i+1))) + ((s_{\text{cb}}(i+2) - s_{\text{cb}}(i)) + (s_{\text{cb}}(i+3) - s_{\text{cb}}(i+1)))^2\right)}{(s_{\text{cb}}(i) - s_{\text{cb}}(i+2))^2 \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i+3))^2} = \frac{120}{36}
$$

$$
\beta_{\text{coef\_}\text{${XYZ}$}}(i+1, 1, 0) \\
= \frac{4 \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i+1))^2 \cdot \left(10 \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i))^2 + (s_{\text{cb}}(i) - s_{\text{cb}}(i-1))^2 + (s_{\text{cb}}(i) - s_{\text{cb}}(i-1)) \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i))\right)}{\left((s_{\text{cb}}(i-1) - s_{\text{cb}}(i+2))^2 \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i+2))^2\right)} = \frac{48}{36}
$$

$$
\beta_{\text{coef\_}\text{${XYZ}$}}(i+1, 1, 1) \\
= \frac{4 \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i+1))^2 \cdot \left((s_{\text{cb}}(i) - s_{\text{cb}}(i+1)) \cdot ((s_{\text{cb}}(i) - s_{\text{cb}}(i-1)) + 20 \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i))) + (2 \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i-1)) + (s_{\text{cb}}(i+1) - s_{\text{cb}}(i))) \cdot (s_{\text{cb}}(i+2) - s_{\text{cb}}(i))\right)}{\left((s_{\text{cb}}(i+1) - s_{\text{cb}}(i-1)) \cdot (s_{\text{cb}}(i-1) - s_{\text{cb}}(i+2))^2 \cdot (s_{\text{cb}}(i+2) - s_{\text{cb}}(i))\right)} = -\frac{60}{36}
$$

$$
\beta_{\text{coef\_}\text{${XYZ}$}}(i+1, 1, 2) \\
= \frac{4 \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i+1))^2 \cdot \left(10 \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i))^2 + (s_{\text{cb}}(i+1) - s_{\text{cb}}(i)) \cdot (s_{\text{cb}}(i+2) - s_{\text{cb}}(i+1)) + (s_{\text{cb}}(i+2) - s_{\text{cb}}(i+1))^2\right)}{\left((s_{\text{cb}}(i-1) - s_{\text{cb}}(i+1))^2 \cdot (s_{\text{cb}}(i-1) - s_{\text{cb}}(i+2))^2\right)} = \frac{48}{36}
$$

$$
\beta_{\text{coef\_}\text{${XYZ}$}}(i+1, 2, 0) \\
= \frac{4 \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i+1))^2 \cdot \left(12 \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i))^2 + \left((s_{\text{cb}}(i) - s_{\text{cb}}(i-2)) + (s_{\text{cb}}(i) - s_{\text{cb}}(i-1))\right)^2 + 3 \cdot \left((s_{\text{cb}}(i) - s_{\text{cb}}(i-2)) + (s_{\text{cb}}(i) - s_{\text{cb}}(i-1))\right) \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i))\right)}{\left((s_{\text{cb}}(i-2) - s_{\text{cb}}(i+1))^2 \cdot (s_{\text{cb}}(i-1) - s_{\text{cb}}(i+1))^2\right)} = \frac{120}{36}
$$

$$
\beta_{\text{coef\_}\text{${XYZ}$}}(i+1, 2, 1) \\
= \frac{4 \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i+1))^2 \cdot \left(19 \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i))^2 + \left((s_{\text{cb}}(i) - s_{\text{cb}}(i-2)) \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i+1))\right) + 2 \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i-1)) \cdot \left((s_{\text{cb}}(i) - s_{\text{cb}}(i-2)) + (s_{\text{cb}}(i+1) - s_{\text{cb}}(i-1))\right)\right)}{\left((s_{\text{cb}}(i-2) - s_{\text{cb}}(i)) \cdot (s_{\text{cb}}(i-2) - s_{\text{cb}}(i+1))^2 \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i-1))\right)} = -\frac{132}{36}
$$

$$
\beta_{\text{coef\_}\text{${XYZ}$}}(i+1, 2, 2) \\
= \frac{4 \cdot (s_{\text{cb}}(i) - s_{\text{cb}}(i+1))^2 \cdot \left(10 \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i))^2 + (s_{\text{cb}}(i) - s_{\text{cb}}(i-1))^2 + (s_{\text{cb}}(i) - s_{\text{cb}}(i-1)) \cdot (s_{\text{cb}}(i+1) - s_{\text{cb}}(i))\right)}{\left((s_{\text{cb}}(i-2) - s_{\text{cb}}(i))^2 \cdot (s_{\text{cb}}(i-2) - s_{\text{cb}}(i+1))^2\right)} = \frac{48}{36}
$$

## s_weno

Variables:

- $\overline{q}$: v_rs_ws_${XYZ}$

Outputs:

- vL_rs_vf: reconstructed left boundary of cell $i$, corresponding to $q_{i-\frac{1}{2}}^R$
- vR_rs_vf: reconstructed right boundary of cell $i$, corresponding to $q_{i+\frac{1}{2}}^L$

- weno_order = 1
  - no reconstruction
- weno_order = 3

    $$
    \Delta_k = \overline{q}_{i+k+1} - \overline{q}_{i+k} \\
    f_0 = \overline{q}_i+c_{0,0}\Delta_0 = \overline{q}_i - \frac{1}{2}\Delta_0 \\
    f_1 = \overline{q}_i + c_{1,0}\Delta_{-1} = \overline{q}_i - \frac{1}{2}\Delta_{-1} \\
    \beta_0 = a_{0,0}\Delta_0^2 + \epsilon_{weno} = \Delta_0^2 + \epsilon_{weno} \\
    \beta_1 = a_{1,0}\Delta_{-1}^2 + \epsilon_{weno} = \Delta_{-1}^2 + \epsilon_{weno} \\
    \alpha_0 = \frac{d_{0,L}}{\beta_0^2} = \frac{1}{3\beta_0^2}, \alpha_1 = \frac{d_{1,L}}{\beta_1^2} = \frac{2}{3\beta_1^2} \\
    \omega_k = \frac{\alpha_k}{\sum\alpha_k}
    $$

    if mapped_weno

    $$
    \alpha_k = (d_{k,L}(1+d_{k,L}-3\omega_k) + \omega_k^2)\cdot(\frac{\omega_k}{d_{k,L}^2 + \omega_k(1-2d_{k,L})}) \\
    \omega_k = \frac{\alpha_k}{\sum\alpha_k}
    $$

    vL_rs_vf $= q_{i-\frac{1}{2}}^R= \omega_0f_0 + \omega_1f_1$

    $$
    f_0 = \overline{q}_i + c_{0,0}^R\Delta_0 = \overline{q}_i + \frac{1}{2}\Delta_0 \\
    f_1 = \overline{q}_i + c_{1,0}^R\Delta_{-1} = \overline{q}_i + \frac{1}{2}\Delta_{-1} \\
    \alpha_0 = d_0^R/\beta_0^2 = \frac{2}{3\beta_0^2}, \alpha_1 = \frac{d_1^R}{\beta_1^2} = \frac{1}{3\beta_1^2} \\
    \omega_k = \frac{\alpha_k}{\sum\alpha_k}
    $$

    if mapped_weno

    $$
    \alpha_k = (d_k^R(1+d_k^R-3\omega_k)+\omega_k^2)\cdot(\frac{\omega_k}{d_{k,R}^2 + \omega_k(1-2d_{k,R})}) \\
    \omega_k = \frac{\alpha_k}{\sum\alpha_k}
    $$

    vR_rs_vf $=q_{i+\frac{1}{2}}^L = \omega_0f_0 + \omega_1f_1$

- weno_order = 5

    $$
    \Delta_k = \overline{q}_{i+k+1} - \overline{q}_{i+k} \\
    f_0 = \overline{q}_i + c_{0,0}\Delta_1 + c_{0,1}\Delta_0 = \overline{q}_i + \frac{2}{6}\Delta_1-\frac{5}{6}\Delta_0 \\
    f_1 = \overline{q}_i + c_{1,0}\Delta_0 + c_{1,1}\Delta_{-1} = \overline{q}_i - \frac{1}{6}\Delta_0 -\frac{2}{6}\Delta_{-1} \\
    f_2 = \overline{q}_i + c_{2,0}\Delta_{-1} + c_{2,1}\Delta_{-2} = \overline{q}_i - \frac{4}{6}\Delta_{-1} + \frac{1}{6}\Delta_{-2} \\
    \beta_0 = a_{0,0}\Delta_1^2  + a_{0,1}\Delta_1\Delta_0 + a_{0,2}\Delta_0^2 + \epsilon_{weno} \\
    = \frac{48}{60}\Delta_1^2 - \frac{132}{36}\Delta_1\Delta_0 + \frac{120}{36}\Delta_0^2 + \epsilon_{weno} \\
    \beta_1 = a_{1,0}\Delta_0^2 + a_{1,1}\Delta_0\Delta_{-1} + a_{1,2}\Delta_{-1}^2 + \epsilon_{weno} \\
    = \frac{48}{36}\Delta_0^2 - \frac{60}{36}\Delta_0\Delta_{-1} + \frac{48}{36}\Delta_{-1}^2 + \epsilon_{weno} \\
    \beta_2 = a_{2,0}\Delta_{-1}^2 + a_{2,1}\Delta_{-1}\Delta_{-2} + a_{2,2}\Delta_{-2}^2 + \epsilon_{weno} \\
    = \frac{120}{36}\Delta_{-1}^2 - \frac{132}{36}\Delta_{-1}\Delta_{-2} + \frac{48}{36}\Delta_{-2}^2 + \epsilon_{weno} \\
    \alpha_0 = \frac{d_{0,L}}{\beta_0^2} = \frac{2}{20\beta_0^2} \\
    \alpha_1 = \frac{d_{1,L}}{\beta_1^2} = \frac{12}{20\beta_1^2} \\
    \alpha_2 = \frac{d_{2,L}}{\beta_2^2} = \frac{6}{20\beta_2^2} \\
    \omega_k = \frac{\alpha_k}{\sum\alpha_k}
    $$

    if mapped_weno

    $$
    \alpha_k = (d_{k,L}(1+d_{k,L}-3\omega_k)+\omega_k^2)\cdot(\frac{\omega_k}{d_{k,L}^2  + \omega_k(1-2d_{k,L})}) \\
    \omega_k = \frac{\alpha_k}{\sum\alpha_k}
    $$

    vL_rs_vf $=q_{i-\frac{1}{2}}^R=\sum\omega_kf_k$

    $$
    f_0 = \overline{q}_i + c_{0,0}\Delta_1 + c_{0,1}\Delta_0 = \overline{q}_i - \frac{1}{6}\Delta_1 + \frac{4}{6}\Delta_0 \\
    f_1 = \overline{q}_i + c_{1,0}\Delta_0 + c_{1,1}\Delta_{-1} = \overline{q}_i + \frac{2}{6}\Delta_0 + \frac{1}{6}\Delta_{-1} \\
    f_2 = \overline{q}_i + c_{2,0}\Delta_{-1} + c_{2,1}\Delta_{-2} = \overline{q}_i + \frac{5}{6}\Delta_{-1} - \frac{2}{6}\Delta_{-2} \\
    \alpha_0 = \frac{d_{0,R}}{\beta_0^2} = \frac{6}{20\beta_0^2} \\
    \alpha_1 = \frac{d_{1,R}}{\beta_1^2} = \frac{12}{20\beta_1^2} \\
    \alpha_2 = \frac{d_{2,R}}{\beta_2^2} = \frac{2}{20\beta_2^2} \\
    \omega_k = \frac{\alpha_k}{\sum\alpha_k}
    $$

    if mapped_weno

    $$
    \alpha_k = (d_{k,R}(1+d_{k,R}-3\omega_k)+\omega_k^2)\cdot(\frac{\omega_k}{d_{k,R}^2 + \omega_k(1-2d_{k,R})}) \\
    \omega_k = \frac{\alpha_k}{\sum\alpha_k}
    $$

    vR_rs_vf $=q_{i+\frac{1}{2}}^L = \sum\omega_kf_k$

if mp_weno

`s_preserve_monotonicity`

## s_initialize_weno

transfer global cell-average variables outside to the local variables.

Note: `s_weno` is utilized in `s_reconstruct_cell_boundary_values` from `m_rhs.fpp`, in which is1,is2,is3 have been reset according to dimensional splitting, specially, is1%beg = is1%beg + weno_polyn, is1%end = is1%end - weno_polyn. **Therefore, in this part, we again $\mp$ weno_polyn to completely get the values from -buff_size to m/n/p+buff_size.**

- weno_dir = 1

v_rs_ws_x(k, l, q, j) = v_vf(j)%sf(k, l, q)

- weno_dir = 2

v_rs_ws_y(k, l, q, j) = v_vf(j)%sf(l, k, q)

- weno_dir = 3

v_rs_ws_z(k, l, q, j) = v_vf(j)%sf(q, l, k)

## s_preserve_monotonicity

## s_finalize_weno_module

deallocate v_rs_ws_x, poly_coef_cbL_x, d_cbL_x, beta_coef_x.
