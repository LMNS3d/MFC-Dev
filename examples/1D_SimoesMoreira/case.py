#!/usr/bin/env python3

import math,json

# Pressure
p01 = 726544.2536731864
p02 = 1.0E+02
# liquid dodecane

# density
rho0wl1 = 480.6864033286308
rho0wl2 = 479.8150064108966
# pi infty
piwl = 4.0E+08
# qv
qvwl = -775.269E+03
# qv'
qvpwl = 0.0E0
# cv
cvwl = 1077.7
# gamma
gamwl = 2.35E0
# cp
cpwl = cvwl*gamwl

alpha0wl1 = 1.0E0 - 1.E-08
alpha0wl2 = 1.E-08

# Vapor dodecane

# density
rho0wv1 = 25.9238015527443
rho0wv2 = 0.003568096701843206
# pi infinity
piwv = 0
# qv
qvwv = -237.547E+03
# qv'
qvpwv = -24400.0E0
# cv
cvwv = 1956.45
# gamma
gamwv = 1.025E0
# cp
cpwv = cvwv*gamwv

alpha0wv1 = 1.0E0 - alpha0wl1
alpha0wv2 = 1.0E0 - alpha0wl2

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : -1.E+00,
    'x_domain%end'                 : 1.E+00,
    'm'                            : 1999,
    'n'                            : 0,
    'p'                            : 0,
    'dt'                           : 5.E-07,
    't_step_start'                 : 0,
    't_step_stop'                  : 5000,
    't_step_save'                  : 100,
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 2,
    'model_eqns'                   : 3,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 2,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'T',
    'mixture_err'                  : 'T',
    'relax'        : 'T',      
    'relax_model'  : 6,	       
    'palpha_eps'   : 1.0E-2,   
    'ptgalpha_eps' : 1.0E-2,   
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'weno_Re_flux'                 : 'F',  
    'weno_avg'                     : 'F',
    'mapped_weno'                  : 'F',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -3,
    'bc_x%end'                     : -3,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'F',
    # ==========================================================================
                                            
    # Patch 1: Left state ======================================================
    'patch_icpp(1)%geometry'       : 1,
    'patch_icpp(1)%x_centroid'     : -1.0E+00/2.0E+00,
    'patch_icpp(1)%length_x'       : 1.0E+00,
    'patch_icpp(1)%vel(1)'         : 0.E+00,
    'patch_icpp(1)%pres'           : p01,
    'patch_icpp(1)%alpha_rho(1)'   : rho0wl1*alpha0wl1,
    'patch_icpp(1)%alpha_rho(2)'   : rho0wv1*alpha0wv1,
    'patch_icpp(1)%alpha(1)'       : alpha0wl1,
    'patch_icpp(1)%alpha(2)'       : alpha0wv1,
    # ==========================================================================

    # Patch 2: Right state =====================================================
    'patch_icpp(2)%geometry'       : 1,
    'patch_icpp(2)%x_centroid'     : 1.0E+00/2.0E+00,
    'patch_icpp(2)%length_x'       : 1.0E+00,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    'patch_icpp(2)%vel(1)'         : 0.E+00,
    'patch_icpp(2)%pres'           : p02,
    'patch_icpp(2)%alpha_rho(1)'   : rho0wl2*alpha0wl2,
    'patch_icpp(2)%alpha_rho(2)'   : rho0wv2*alpha0wv2,
    'patch_icpp(2)%alpha(1)'       : alpha0wl2,
    'patch_icpp(2)%alpha(2)'       : alpha0wv2,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.0E+00 / ( gamwl - 1 ),       
    'fluid_pp(1)%pi_inf'           : gamwl * piwl / ( gamwl - 1 ),  
    'fluid_pp(1)%cv'          	   : cvwl,          
    'fluid_pp(1)%qv'        	   : qvwl,	
    'fluid_pp(1)%qvp'          	   : qvpwl,         
    'fluid_pp(2)%gamma'            : 1.0E+00 / ( gamwv - 1 ),       
    'fluid_pp(2)%pi_inf'           : gamwv * piwv / ( gamwv - 1 ),  
    'fluid_pp(2)%cv'          	   : cvwv,          
    'fluid_pp(2)%qv'        	   : qvwv,  	
    'fluid_pp(2)%qvp'          	   : qvpwv,	
    # ==========================================================================
}))

# ==============================================================================
