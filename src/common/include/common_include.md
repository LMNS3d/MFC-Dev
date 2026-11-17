# inline_conversion.fpp (sound speed calculation)

## calculate c^2 for different model

- alt_soundspeed (Kapila's five-equation model)

$$
gammas_k = \Gamma_k = \frac{1}{\gamma_k - 1} \\
pi_{infs,k} = \Pi_k = \frac{\gamma_k p_k^{\infty}}{\gamma_k - 1} \\
blkmod_k = \frac{(\Gamma_k+1)p + \Pi_k}{\Gamma_k} = \gamma_k(p+p_k^{\infty}) \\
c^2 = \frac{1}{\rho(\frac{\alpha_1}{blkmod_1}+\frac{\alpha_2}{blkmod_2})} = \frac{1}{\rho(\frac{\alpha_1}{\gamma_1(p+p_1^{\infty})}+\frac{\alpha_2}{\gamma_2(p+p_2^{\infty})})}
$$

- model_eqns == 3 (Saurel's six-equation model)

$$
c^2 = \frac{\sum{\alpha_k(1/\Gamma_k+1)(p+\Pi_k/(\Gamma_k+1))}}{\rho} = \frac{\sum{\alpha_k\gamma_k(p+p_k^{\infty})}}{\rho}
$$

- Allaire's five-equation model

$$
gamma = \Gamma = \frac{1}{\gamma - 1} = \sum{\alpha_k\Gamma_k} \\
pi_{inf} = \Pi = \frac{\gamma p^{\infty}}{\gamma - 1} = \sum{\alpha_k\Pi_k} \\
H = (E+p)/\rho = \frac{(\Gamma+1) p + \Pi +\frac{1}{2}\rho|u|^2}{\rho} \\
c^2 = \frac{H-\frac{1}{2}|u|^2}{\Gamma} = \frac{(\Gamma+1) p + \Pi}{\Gamma\rho} = \frac{\gamma(p+p^{\infty})}{\rho}
$$

## calculate c with modification if required

- mixture_err 
  - if c<0, set c to be a small positive number.

# macros.fpp

Use preprocessor macros to define operations related to setting up variables for use with OpenACC.

