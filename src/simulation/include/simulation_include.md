# inline_riemann.fpp

## arithmetic_avg

$$
\rho_{avg} = \frac{\rho_L + \rho_R}{2} \\
u_{avg}^2 = (\frac{u_L + u_R}{2})^2 + (\frac{v_L + v_R}{2})^2 \\
H_{avg} = \frac{H_L + H_R}{2} \\
\Gamma_{avg} = \frac{\Gamma_L + \Gamma_R}{2}
$$

## roe_avg

$$
\rho_{avg} = \sqrt{\rho_L\rho_R} \\
u_{avg}^2 = \frac{(\sqrt{\rho_L}u_L + \sqrt{\rho_R}u_R)^2}{(\sqrt{\rho_L} + \sqrt{\rho_R})^2} + \frac{(\sqrt{\rho_L}v_L + \sqrt{\rho_R}v_R)^2}{(\sqrt{\rho_L} + \sqrt{\rho_R})^2} \\
H_{avg} = \frac{\sqrt{\rho_L}H_L + \sqrt{\rho_R}H_R}{\sqrt{\rho_L} + \sqrt{\rho_R}} \\
\Gamma_{avg} = \frac{\sqrt{\rho_L}\Gamma_L + \sqrt{\rho_R}\Gamma_R}{\sqrt{\rho_L} + \sqrt{\rho_R}} \\
\rho_{avg} = \sqrt{\rho_L\rho_R} \\
u_{avg}^2 = \frac{(\sqrt{\rho_L}u_L + \sqrt{\rho_R}u_R)^2}{(\sqrt{\rho_L} + \sqrt{\rho_R})^2}
$$

Note that this part seems to be redundant, especially $u_{avg}^2$, which is first calculated by adding all the directions' velocity, but then recalculated for only the split direction.

## compute_average_state

- avg_state = 1: `@:roe_avg()`
- avg_state = 2: `@:arithmetic_avg()`

# case.fpp

nothing inside, utilized for GPU-oriented optimization.

