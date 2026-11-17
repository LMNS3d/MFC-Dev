# Hardcode for initial conditions (physical parameters)

This part is mainly for those cases which need special initial physical parameters setting in each region. For example, the case where the initial velocity is related to location
$$
\psi(i,j) = \sin^2(\pi x)\sin^2(\pi y)\cos(\pi t/2)/\pi \\
u = -\frac{\partial \psi}{\partial y}, v = \frac{\partial \psi}{\partial x}
$$

## 1dHardcodedIC.fpp

...

## 2dHardcodedIC.fpp

1. declare intermediate variables utilized during initial conditons setting

2. select self-defined conditions

   - hcid: index of the defined setting

- hcid = 200

    if $y \le (-x^3 + 1)^{1/3}$
    $$
    \alpha_1 = 10^{-8} \\
    \alpha_2 = 1 - 10^{-8} \\
    \alpha_1\rho_1 = 10^{-8} \times 1000 \\
    \alpha_2\rho_2 = (1 - 10^{-8})\times 1 \\
    p = 1000
    $$

- hcid = 202 (Gresho vortex)

    $$
    r = \sqrt{(x- 0.5)^2 + (y-0.5)^2} \\
    r_{max} = 0.2 \\
    \gamma_1 = 1 + 1/\Gamma_1 \\
    u_{max} = 2\pi r_{max}v \\
    p_0 = u_{max}^2(\frac{1}{\gamma_1 v^2} - 0.5)
    $$

    if $r< r_{max}$
    $$
    u = -(y- 0.5)u_{max}/r_{max} \\
    v = (x - 0.5)u_{max}/r_{max} \\
    p = p_0 + u_{max}^2\frac{(r/r_{max})^2}{2} \\
    $$

    elseif $r_{max} \le r < 2r_{max}$
    $$
    u = -(y - 0.5)u_{max}(2 - r/r_{max})/r \\
    v = (x - 0.5)u_{max}(2 - r/r_{max})/r \\
    p = p_0 + u_{max}^2(\frac{(r/r_{max})^2}{2} + 4(1-r/r_{max} + \ln({r/r_{max}}))) \\
    $$

    else $r\ge 2r_{max}$
    $$
    u = 0 \\
    v = 0 \\
    p = p_0 + u_{max}^2(-2 + 4\ln2)
    $$

- hcid = 203 (similar to 202 but with density correction)

    1. same procedure as hcid = 202
    2. $\alpha_1\rho_1 = p^{\frac{1}{\gamma_1}}$

## 3dHardcodedIC.fpp

...

