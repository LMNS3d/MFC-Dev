# m_phase_change.fpp

This module is used to relax the 6-eqn model towards pressure and temperature, and Gibbs free energies equilibrum through an infinitely fast procedure.
In general, the pressure equilibrum well be maintained at any place, while temperature and Gibbs energy equilibrium are maintained at the material interface.

The procedure is performed ```after``` RK3 time integration.

- input: q_cons_vf
  - before consider phase change, the two phases have different temperature and different Gibbs energies.
- output: q_cons_vf
  - after this procedure, the two phases will share the same temperature and may have the same Gibbs energy.

parameters:

- $p_{Cr}$ = 4.97e7, $T_{Cr}$ = 385.05 + 273.15
  - Out of this region, water will have very strange state, which is not liquid and vapor, and the model is not suitable anymore.
- mixM = 1e-8
  - threshold for 'mixture cell', if Y < mixM, we consider it to be single phase, therefore the phase change is not activated
- lp = 1 index of liquid phase
- vp = 2 index of vapor phase

## s_initialize_phasechange_module (calculate parameter for g1=g2 function)

$$
gs\_min_k = \gamma_k = 1/\Gamma_k +1 \\
ps\_inf_k = p^{\infty}_k = \Pi_k/(1+\Gamma_k) \\
C_{p,k} = \gamma_kC_{v,k}
$$
$$
A = \frac{\gamma_1 C_{v,1} - \gamma_2 C_{v,2} + q^{\prime}_2 - q^{\prime}_1}{C_{v,2}(\gamma_2 - 1)} \\
B = \frac{q_1 - q_2}{C_{v,2}(\gamma_2 - 1)} \\
C = \frac{\gamma_2C_{v,2} - \gamma_1C_{v,1}}{C_{v,2}(\gamma_2 - 1)}\\
D = \frac{C_{v,1}(\gamma_1 - 1)}{C_{v,2}(\gamma_2 - 1)}
$$

## s_infinite_relaxation_k (main procedure for relaxation)

Overheated vapor (OV) and subcooled liquid (SL) are stable states, and phase change does not happen within these states.

variables:

- $p_S,T_S$: equilibrium pressure and temperature for mixture
- $p_{SOV},T_{SOV}$: overheated vapor's pressure and temperature
- $p_{SSL},T_{SSL}$: subcooled liquid's pressure and temperature
- $T_{SatOV},T_{SatSL}$: Saturation Temperatures for OV and SL (under given pressure)
- MCT: partial density correction parameter
- TvF: total volume fraction

1. calculate some unchanged variables

    $$
    \rho = \sum{\alpha_k\rho_k} \\
    TvF = \sum{\alpha_k} \\
    rM = \alpha_1\rho_1 +\alpha_2\rho_2
    $$
    - call ```s_correct_partial_densities```
    $$
    m1 = \alpha_1\rho_1 \\
    m2 = \alpha_2\rho_2 \\
    dynE = \frac{1}{2}\rho |u|^2 \\
    \rho e = E - \frac{1}{2}\rho |u|^2
    $$
    Note that $\rho, \rho e, rM$ will not change throughout phase-change process.

2. pT-equilibrium (thermal relaxation T1=T2)

    - call ```s_infinite_pt_relaxation_k```

    **Note that the thermal relaxation is not restricted within mixture region in the code, at least I do not see any procedure related to it.** $\alpha_k\rho_k$ will not change during this step.

3. pTg-equilibrium if needed

    During this step, mass transition happen between liquid and vapor phase, therefore $\alpha_k\rho_k$ will change during this step.

    general conditions
    - relax_model == 6
    - $Y_{1,rM} = \frac{\alpha_1\rho_1}{rM}> 10^{-8}$, $Y_{2,rM} = \frac{\alpha_2\rho_2}{rM}> 10^{-8}$
      - make sure the cell is within mixture region
    - $p_S<p_{Cr}, T_S <T_{Cr}$
      - make sure the physical conditions are located in reasonable phase-diagram regions

    assume final solution to be OV
    $$
    \alpha_1\rho_1 = 10^{-8}rM \\
    \alpha_2\rho_2 = (1-10^{-8})rM
    $$
    - call ```s_infinite_pt_relaxation_k``` to get $p_{SOV},T_{SOV}$
    - call ```s_TSat``` to get $T_{SatOV}$
    - if $T_{SOV} > T_{SatOV}$, accept this assumption: $T_S = T_{SOV},p_S = p_{SOV}$

    assume final solution to be SL
    $$
    \alpha_1\rho_1 = (1-10^{-8})rM \\
    \alpha_2\rho_2 = 10^{-8}rM
    $$
    - call ```s_infinite_pt_relaxation_k``` to get $p_{SSL},T_{SSL}$
    - call ```s_TSat``` to get $T_{SatSL}$
    - if $T_{SSL} < T_{SatSL}$, accept this assumption: $T_S = T_{SSL},p_S = p_{SSL}$

    Otherwise, reject these assumptions
    $$
    \alpha_1\rho_1 = m1 \\
    \alpha_2\rho_2 = m2
    $$
    - call ```s_infinite_ptg_relaxation_k``` to calculate $p_S,T_S$

4. calculations after all the relaxation processes

$$
s_k = C_{v,k}\ln{\frac{T_S^{\gamma_k}}{(p_S + p^{\infty}_k)^{\gamma_k-1}}} + q^{\prime}_k \\
h_k = \gamma_kC_{v,k}T_S +q_k \\
g_k = h_k - T_S s_k \\
\rho_k = \frac{p_S + p^{\infty}_k}{(\gamma_k - 1)C_{v,k}T_S} \\
e_k = \frac{(p_S + \gamma_k p^{\infty}_k)}{p_S + p^{\infty}_k} C_{v,k}T_S + q_k \\
\alpha_k = \frac{\alpha_k\rho_k}{\rho_k} \\
\alpha_k\rho_k e_k = (\alpha_k\rho_k)e_k \\
\rho s = \sum{(\alpha_k\rho_k) s_k}
$$

## s_infinite_pt_relaxation_k (thermal relaxation)

After this relaxation, the two phases will have the same temperature $T_S$ and of course the same pressure $p_S$. Again, **the relaxation step is not restricted within interface region by this code.**

Input:

- ```q_cons_vf```

Output:

- $p_S$
- $T_S$

Note that ```q_cons_vf``` is not changed during this subroutine. Therefore, it should be updated outside.

Variables:

- $p\_infpT_k = p^{\infty}_k$: stiffness for the participating fluids under pT-equilibrium
- MFL: flag to tell wheter the fluid is pure gas (0), pure liquid (1), or a mixture (2)
- $gp=g(p_S), gpp=g^{\prime}(p_S)$
- $hp=h(p_S)$ is an auxiliary function for better convergence
- $p_O$: initial pressure guess and old pressure during calculation
- $palpha\_eps=\epsilon_p$: error parameter for calculating pressure by Newton's method

$$
mCP = \rho C_p = \sum{(\alpha_k\rho_k)\gamma_kC_{v,k}} \\
mQ = \rho q = \sum{(\alpha_k\rho_k) q_k}
$$

1. check energy constraint
   - $\rho e -\rho q < \min(p^{\infty}_k)$, where $\min(p^{\infty}_k)$ usually equals 0
   - if MFL = 0 or 1, which means that the fluid is pure phase, we set $p_S = 0,T_S=0$ and finish the subroutine

2. solve the ODE of the relaxation process by Newton's method

    1. initially set $p_O = 0,p_S = 10^4$

    2. iteration to calculate $p_S$

    Note: This Newton's method is based on the equation $1 = T_S\sum\frac{\alpha_k}{T_k}$, since the temperature is the same after relaxation.

    - $|p_S - p_O|>\epsilon_p$
    - $|p_S - p_O|/p_O >\epsilon_p/10^4$
    $$
    p_O = p_S \\
    g(p_S) = \sum{\frac{(\gamma_k-1)C_{v,k}\alpha_k\rho_k(\rho e -\rho q + p_S)}{\rho C_p(p_S+p^{\infty}_k)}} \\
    g^{\prime}(p_S) = \sum{\frac{(\gamma_k-1)C_{v,k}\alpha_k\rho_k(p^{\infty}_k-\rho e +\rho q)}{\rho C_p(p_S + p^{\infty}_k)^2}} \\
    h(p_S) = \frac{1}{\rho e + p_S - \rho q} + \frac{1}{p_S + \min(p^{\infty}_k)} \\
    p_S = p_O + \frac{(1 - g(p_S)) / g^{\prime}(p_S)}{1 - \left( \frac{1 - g(p_S) + |1 - g(p_S)|}{2 \cdot g^{\prime}(p_S)} \right) \cdot h(p_S)}
    $$

    3. compute $T_S$ based on $p_S$
    $$
    T_S = \frac{\rho e - \rho q +p_S}{\rho C_p}
    $$

## s_infinite_ptg_relaxation (phase change around interface g1=g2)

Variables:

- $p\_infpTg_k = p^{\infty}_k$
- $O_m=10^{-3}$: underrelaxation factor 
- $mCP = \rho C_p$
- $mQ=\rho q$
- $mCPD = \alpha_3\rho_3 C_{p,3}$: for non-reacting phase (only calculated when considering three-phase flow)
- $mQD=\alpha_3\rho_3 q_3$: for non-reacting phase (similar to variable above)
- $R_{2D}$: residuals for two equations
- $ptgalpha\_eps = \epsilon_{ptg}$: error parameter for the allowable residuals
- $DeltamP = \Delta X$: corrections for two unknowns
- $Jac = J, InvJac = J^{-1}, TJac = J^T$

1. improve initial condition of pressure if needed
    - $p_S <0$ and $\alpha_1\rho_1 +\alpha_2\rho_2 >\frac{\rho e - \frac{\gamma_1p^{\infty}_1}{\gamma_1 - 1}}{q_1}$
  
        or

    - $0\leq p_S<10^{-1}$

    $p_S = 10^4$ if the above condition satisfied. **This part not understand.**

2. solve ODE of thermal-chemical relaxation by Newton's method (two unknown variables $\alpha_1\rho_1,p_S$ given two equations $F_1,F_2$).

    1. loop condition
       - $\sqrt{R_{2D,1}^2 + R_{2D,2}^2} > \epsilon_{ptg}$
       - $\frac{\sqrt{R_{2D,1}^2 + R_{2D,2}^2}}{\rho e}>\epsilon_{ptg}/10^6$

    2. calculation within the loop
        $$
        \rho C_p = \sum{(\alpha_k\rho_k)\gamma_k C_{v,k}} \\
        \rho q = \sum{(\alpha_k\rho_k)q_k} \\
        mCPD = (\alpha_3\rho_3)\gamma_3 C_{v,3} \\
        mQD = (\alpha_3\rho_3) q_3 \\
        mCVGP = \frac{\alpha_3\rho_3(\gamma_3 - 1)C_{v,3}}{p_S + p^{\infty}_3} \\
        mCVGP2 = \frac{\alpha_3\rho_3(\gamma_3 - 1)C_{v,3}}{(p_S + p^{\infty}_3)^2}
        $$

        call ```s_compute_jacobian_matrix``` to get $J,J^{-1},J^{T}$
        $$
        \Delta X = - J^{-1}R_{2D} \\
        \alpha_1\rho_1 = \alpha_1\rho_1 + O_m\Delta X_1 \\
        \alpha_2\rho_2 = \alpha_2\rho_2 - O_m\Delta X_1 \\
        p_S = p_S + O_m\Delta X_2
        $$

        call ```s_compute_pTg_residue``` to calculate $R_{2D}$

    3. calculate $T_S$ based on $p_S$
        $$
        T_S = \frac{\rho e - \rho q +p_S}{\rho C_p}
        $$

## s_correct_partial_densities (liquid and vapor phases correction)

parameter:

- $MCT = 2 \times 10^{-8}$: partial density correction parameter
- $rM = \alpha_1\rho_1 + \alpha_2\rho_2$: sum of reacting masses

1. correction for negative total reacting masses
    - $rM < 0$
    - $\alpha_1\rho_1 \geq -10^{-8}$
    - $\alpha_2\rho_2 \geq -10^{-8}$
    $$
    \alpha_1\rho_1 = 0 \\
    \alpha_2\rho_2 = 0 \\
    rM = 0
    $$

2. correction for partial densities of reacting fluids
    - if $\alpha_1\rho_1 < 0$
    $$
    \alpha_1\rho_1 = MCT \cdot rM \\
    \alpha_2\rho_2 = (1-MCT) \cdot rM
    $$
    - if $\alpha_2\rho_2 < 0$
    $$
    \alpha_1\rho_1 = (1-MCT) \cdot rM \\
    \alpha_2\rho_2 = MCT \cdot rM
    $$

## s_compute_jacobian_matrix (used in pTg-equilibrium)

$$
F = (F_1,F_2)^T \\
X = (m_l,p_S)^T \\
J = \frac{dF}{dX} = \begin{pmatrix}
\frac{\partial F_1}{\partial m_l} & \frac{\partial F_1}{\partial p_S} \\
\frac{\partial F_2}{\partial m_l} & \frac{\partial F_2}{\partial p_S}
\end{pmatrix}
$$

where

- $F_1 = F_1(m_l,p_S)=0$: Gibbs free energy quality
- $F_2 = F_2(m_l,p_S)=0$: enforcement of thermodynamic energy $\rho e$
- $m_l = \alpha_1\rho_1$: mass of reacting liquid (changed)
- $m_T = \alpha_1\rho_1 +\alpha_2\rho_2$: total reacting masses (unchanged)

Use $m_l,p_S$ to describe all other variables or functions
$$
T_S = \frac{1}{m_T\frac{(\gamma_2 - 1)C_{v,2}}{p_S +p^{\infty}_2} + m_l(\frac{(\gamma_1 - 1)C_{v,1}}{p_S + p^{\infty}_1} - \frac{(\gamma_2 - 1)C_{v,2}}{p_S + p^{\infty}_2}) + mCVGP}
$$
$$
\begin{aligned}
dFdT = & \frac{\partial F_1(T_S, p_S)}{\partial T_S} \\
=& - (\gamma_1 C_{v,1} - \gamma_2 C_{v,2}) \ln T_S \\
& - (q^{\prime}_1 - q^{\prime}_2) \\
& + (\gamma_1 - 1) C_{v,1} \ln (p_S + p^{\infty}_1) \\
& - (\gamma_2 - 1) C_{v,2} \ln (p_S + p^{\infty}_2)
\end{aligned}
$$
$$
\begin{aligned}
dTdm = & \frac{\partial T_S(m_l,p_S)}{\partial m_l} \\
=& - \left( \frac{C_{v,1} (\gamma_1 - 1)}{p_S + p^{\infty}_1} - \frac{C_{v,2} (\gamma_2 - 1)}{p_S + p^{\infty}_2} \right) T_S^2
\end{aligned}
$$
$$
\begin{aligned}
dTdp = & \frac{\partial T_S(m_l,p_S)}{\partial p_S} \\
=& ( m_T \frac{C_{v,2} (\gamma_2 - 1)}{(p_S + p^{\infty}_2)^2} \\
&+ m_l ( \frac{C_{v,1} (\gamma_1 - 1)}{(p_S + p^{\infty}_1)^2} - \frac{C_{v,2} (\gamma_2 - 1)}{(p_S + p^{\infty}_2)^2}) \\
&+ mCVGP2 ) T_S^2
\end{aligned}
$$

$$
\begin{aligned}
J_{1,1} = &\frac{\partial F_1(m_l,p_S)}{\partial m_l} \\ 
=& \frac{\partial F_1(T_S, p_S)}{\partial T_S} \cdot \frac{\partial T_S(m_l,p_S)}{\partial m_l} \\
=& dFdT \cdot dTdm
\end{aligned} \\
\begin{aligned}
J_{1,2} = &\frac{\partial F_1(m_l,p_S)}{\partial p_S} \\
=& \frac{\partial F_1(T_S, p_S)}{\partial T_S} \cdot \frac{\partial T_S(m_l,p_S)}{\partial p_S} + \frac{\partial F_1(T_S, p_S)}{\partial p_S} \\
= & dFdT \cdot dTdp + T_S\left(\frac{(\gamma_1 - 1)C_{v,1}}{p_S + p^{\infty}_1} - \frac{(\gamma_2 - 1)C_{v,2}}{p_S + p^{\infty}_2}\right)
\end{aligned}
$$

$$
\begin{aligned}
J_{2,1} = & \frac{\partial F_2(m_l,p_S)}{\partial m_l} \\
= & q^{\prime}_2 - q^{\prime}_1 \\
& + \frac{\gamma_2 C_{v,2} - \gamma_1 C_{v,1}}{m_l \left( \frac{C_{v,1} (\gamma_1 - 1)}{p_S + p^{\infty}_1} - \frac{C_{v,2} (\gamma_2 - 1)}{p_S + p^{\infty}_2} \right) + m_T \frac{C_{v,2} (\gamma_2 - 1)}{p_S + p^{\infty}_2} + mCVGP} \\
& - \left( m_l (\gamma_2 C_{v,2} - \gamma_1 C_{v,1}) - m_T \gamma_2 C_{v,2} - mCPD \right) \\
& \cdot \frac{ \left( \frac{C_{v,1} (\gamma_1 - 1)}{p_S + p^{\infty}_1} - \frac{C_{v,2} (\gamma_2 - 1)}{p_S + p^{\infty}_2} \right)}{(m_l \left( \frac{C_{v,1} (\gamma_1 - 1)}{p_S + p^{\infty}_1} - \frac{C_{v,2} (\gamma_2 - 1)}{p_S + p^{\infty}_2} \right) + m_T \frac{C_{v,2} (\gamma_2 - 1)}{p_S + p^{\infty}_2} + mCVGP)^2}  \\
=& q^{\prime}_2 - q^{\prime}_1 \\
& +T_S(\gamma_2 C_{v,2} - \gamma_1 C_{v,1}) \\
& - \left( m_l (\gamma_2 C_{v,2} - \gamma_1 C_{v,1}) - m_T \gamma_2 C_{v,2} - mCPD \right) \\
& \cdot T_S^2 \left( \frac{C_{v,1} (\gamma_1 - 1)}{p_S + p^{\infty}_1} - \frac{C_{v,2} (\gamma_2 - 1)}{p_S + p^{\infty}_2} \right) \\
=& q^{\prime}_2 - q^{\prime}_1 +T_S(\gamma_2 C_{v,2} - \gamma_1 C_{v,1}) \\
& + \left( m_l (\gamma_2 C_{v,2} - \gamma_1 C_{v,1}) - m_T \gamma_2 C_{v,2} - mCPD \right) \cdot dTdm
\end{aligned}
$$

$$
\begin{aligned}
J_{2,2} = & \frac{\partial F_2(m_l,p_S)}{\partial p_S} \\
= & 1 + \left( m_l (\gamma_2 C_{v,2} - \gamma_1 C_{v,1}) - m_T \gamma_2 C_{v,2} - mCPD \right) \\
& \cdot \frac{( m_T \frac{C_{v,2} (\gamma_2 - 1)}{(p_S + p^{\infty}_2)^2} + m_l ( \frac{C_{v,1} (\gamma_1 - 1)}{(p_S + p^{\infty}_1)^2} - \frac{C_{v,2} (\gamma_2 - 1)}{(p_S + p^{\infty}_2)^2}) + mCVGP2 )}{(m_l \left( \frac{C_{v,1} (\gamma_1 - 1)}{p_S + p^{\infty}_1} - \frac{C_{v,2} (\gamma_2 - 1)}{p_S + p^{\infty}_2} \right) + m_T \frac{C_{v,2} (\gamma_2 - 1)}{p_S + p^{\infty}_2} + mCVGP)^2} \\
= & 1 + \left( m_l (\gamma_2 C_{v,2} - \gamma_1 C_{v,1}) - m_T \gamma_2 C_{v,2} - mCPD \right)\cdot dTdp
\end{aligned}
$$

$$
J^{-1} = \frac{1}{J_{1,1}J_{2,2}-J_{1,2}J_{2,1}} \begin{pmatrix}
J_{2,2} & -J_{1,2} \\
-J_{2,1} & J_{1,1}
\end{pmatrix}
$$

$$
J^T = \begin{pmatrix}
J_{1,1} & J_{2,1} \\
J_{1,2} & J_{2,2}
\end{pmatrix}
$$

## s_compute_pTg_residue (used in pTg-equilibrium)

Use $m_l = \alpha_1\rho_1,p_S$ to calculate all the variables
$$
m_T = \alpha_1\rho_1 + \alpha_2\rho_2 \\
T_S = \frac{1}{m_T\frac{(\gamma_2 - 1)C_{v,2}}{p_S +p^{\infty}_2} + m_l(\frac{(\gamma_1 - 1)C_{v,1}}{p_S + p^{\infty}_1} - \frac{(\gamma_2 - 1)C_{v,2}}{p_S + p^{\infty}_2}) + mCVGP}
$$

- Gibbs free energy equality condition

$$
\begin{aligned}
    R_{2D,1} =& F_1(m_l,p_S) \\
    = & T_S [(\gamma_1C_{v,1}-\gamma_2C_{v,2})(1-\ln (T_S)) - (q^{\prime}_1-q^{\prime}_2) \\
    & + (\gamma_1 - 1)C_{v,1}\ln (p_S + p^{\infty}_1) \\
    & - (\gamma_2 - 1)C_{v,2}\ln (p_S + p^{\infty}_2) ] + q_1 - q_2
\end{aligned}
$$

- Constant total internal energy process condition

$$
\begin{aligned}
    R_{2D,2} =& F_2(m_l,p_S) \\
    = & \rho e + p_S \\
    & + m_l(q_2-q_1) - m_Tq_2 - mQD \\
    & + \frac{\left( m_l (\gamma_2 C_{v,2} - \gamma_1 C_{v,1}) - m_T \gamma_2 C_{v,2} - mCPD \right)}{m_T\frac{(\gamma_2 - 1)C_{v,2}}{p_S +p^{\infty}_2} + m_l(\frac{(\gamma_1 - 1)C_{v,1}}{p_S + p^{\infty}_1} - \frac{(\gamma_2 - 1)C_{v,2}}{p_S + p^{\infty}_2}) + mCVGP} \\
    = & \rho e +p_S + m_l(q_2-q_1) - m_Tq_2 - mQD \\
    & + \left( m_l (\gamma_2 C_{v,2} - \gamma_1 C_{v,1}) - m_T \gamma_2 C_{v,2} - mCPD \right) T_S
\end{aligned}
$$

## s_TSat (Saturation temperature calculation under given pressure)

- $T_{Sat}$: saturation temperature (output)
- $T_{SIn}$: equilibrium temperature (input)
- $p_{Sat}$: saturation pressure (input)
- $O_m = 10^{-3}$: underrelaxation factor

If $p_{Sat} = 0$ and $T_{SIn} = 0$ (energy constraint not satisfied, which means unphysical)
$$
T_{Sat} = 0
$$

otherwise, use saturature curve to calculate $T_{Sat}$ with given $p_{Sat}$, by Newton's method

1. initial estimation for $T_{Sat}$

$$
T_{Sat} = T_{SIn}
$$

2. iteration over Newton's solver
    - loop condition
    $$
    |F(T_{Sat})| > \epsilon_{ptg}
    $$
    - residual calculation
    $$
    \begin{aligned}
        F(T_{Sat}) =& T_{Sat} [(\gamma_1C_{v,1}-\gamma_2C_{v,2})(1-\ln (T_{Sat})) - (q^{\prime}_1-q^{\prime}_2) \\
        & + (\gamma_1 - 1)C_{v,1}\ln (p_{Sat} + p^{\infty}_1) \\
        & - (\gamma_2 - 1)C_{v,2}\ln (p_{Sat} + p^{\infty}_2) ] + q_1 - q_2
    \end{aligned}
    $$
    - derivative calculation
    $$
    \begin{aligned}
        dFdT = &\frac{dF(T_{Sat})}{dT_{Sat}} = F^{\prime}(T_{Sat}) \\
        =& -(\gamma_1C_{v,1}-\gamma_2C_{v,2})\ln (T_{Sat}) - (q^{\prime}_1-q^{\prime}_2) \\
        & + (\gamma_1 - 1)C_{v,1}\ln (p_{Sat} + p^{\infty}_1) \\
        & - (\gamma_2 - 1)C_{v,2}\ln (p_{Sat} + p^{\infty}_2)
    \end{aligned}
    $$
    - value updating
    $$
    T_{Sat} = T_{Sat} - \frac{F(T_{Sat})}{F^{\prime}(T_{Sat})}\cdot O_m
    $$
