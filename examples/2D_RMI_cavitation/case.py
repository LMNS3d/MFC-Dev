#!/usr/bin/env python3
import math, json
## 1 FOR BACKGROUND, 2 FOR SHOKED STATE, 3 FOR WATER REGION (WHEN NEEDED)
# Pressure
p01 = 1.E+05
p02 = 1.E+05
p03 = 1.E+05

## FLUID PROPERTIES FOR EACH PATCH


# liquid water

# pi infty
piwl = 1.0E+09
# qv
qvwl = -1167000
# qv'
qvpwl = 0.0E0
# cv
cvwl = 1816
# cp
cpwl = 4267
# gamma
gamwl = 2.35E0

# density
rho0wl1 = 1000.E+00
rho0wl2 = 1000.E+00
rho0wl3 = 1000.E+00

# volume fraction
awl1 = 1.0E+00 - 1.0E-06
awl2 = 1.0E-06
awl3 = 1.0E-06

## Vapor water

# pi infinity
piwv = 0
# qv
qvwv = 2030000
# qv'
qvpwv = -23400
# cv
cvwv = 1040
# cp
cpwv = 1487
# gamma
gamwv = 1.43E0

# density
rho0wv1 = 100.E+00
rho0wv2 = 100.E+00
rho0wv3 = 100.E+00

# volume fraction
awv1 = 1.0E-06
awv2 = 1.0E+00 - 1.0E-06
awv3 = 1.0E+00 - 1.0E-06

# discretization parameters for the x direction
# number of elements
Nx      = 599

# domain boundaries
xb = 0.00
xe = 3.00

# typical cell size
dx	= ( xe - xb ) / Nx

# discretization parameters for the y direction
# number of elements
Ny      = 199

# domain boundaries
yb = 0.0
ye = 1.0

# typical cell size
dy	= ( ye - yb ) / Ny



# Configuring case dictionary ==================================================
print(json.dumps({
    # Logistics ================================================
    'run_time_info': 'T',
    # ==========================================================
    # Computational Domain Parameters ==========================
    'x_domain%beg' : xb,        
    'x_domain%end' : xe,        
    'y_domain%beg' : yb,        
    'y_domain%end' : ye,        
    'm'            : Nx,        
    'n'            : Ny,        
    'p'            : 0,         
    'dt'           : 5.E-07,        
    't_step_start' : 0,         
    't_step_stop'  : 20000,        
    't_step_save'  : 500,         
    # ==========================================================
    # Simulation Algorithm Parameters ==========================
    'num_patches'  : 3,        
    'model_eqns'   : 3,     
    'alt_soundspeed'               : 'F',   
    'num_fluids'   : 2,        
    'adv_alphan'   : 'T',      
    'mpp_lim'      : 'T',      
    'mixture_err'  : 'T',      
    'relax'        : 'T',      
    'relax_model'  : 6,        
    'palpha_eps'   : 1.0E-2,   
    'ptgalpha_eps' : 1.0E-2,   
    'time_stepper' : 3,        
    'weno_order'   : 5,        
    'weno_eps'     : 1.0E-16,
    'mapped_weno'  : 'F',      
    'null_weights' : 'F',      
    'mp_weno'      : 'F',      
    'riemann_solver'    : 2,   
    'wave_speeds'  : 1,        
    'avg_state'    : 2,        
    'bc_x%beg'     : -2,       
    'bc_x%end'     : -3,       
    'bc_y%beg'     : -2,       
    'bc_y%end'     : -2,       
    # ==========================================================
    # Formatted Database Files Structure Parameters ============
    'format'       : 1,        
    'precision'    : 2,        
    'prim_vars_wrt':'T',       
    'parallel_io'  :'F',       
    # ==========================================================
    # Patch 1 - Background =====================================
    'patch_icpp(1)%geometry'       : 3,     
    'patch_icpp(1)%x_centroid'     : ( xe + xb ) * 2 / 4,           
    'patch_icpp(1)%y_centroid'     : ( ye + yb ) * 2 / 4,           
    'patch_icpp(1)%length_x'       : ( xe - xb ) * 2 / 2,           
    'patch_icpp(1)%length_y'       : ( ye - yb ) * 2 / 2,           
    'patch_icpp(1)%vel(1)'         : -200.0E+00,   
    'patch_icpp(1)%vel(2)'         : 0.0E+00,   
    'patch_icpp(1)%pres'           : p01,  	
    'patch_icpp(1)%alpha_rho(1)'   : awl1 * rho0wl1,           	
    'patch_icpp(1)%alpha_rho(2)'   : awv1 * rho0wv1,
    'patch_icpp(1)%alpha(1)'       : awl1,   	
    'patch_icpp(1)%alpha(2)'       : awv1,
    # ==========================================================
    # Patch 2 - vapor bubble ==================================
    'patch_icpp(2)%geometry'       : 2,     	
    'patch_icpp(2)%x_centroid'     : 1.2E+00,       
    'patch_icpp(2)%y_centroid'     : ( ye + yb ) * 2 / 4,       
    'patch_icpp(2)%radius'         : 0.6E00,
    'patch_icpp(2)%vel(1)'         : -200.0E00,   	
    'patch_icpp(2)%vel(2)'         : 0.0E00,   	    
    'patch_icpp(2)%pres'           : p02,    	
    'patch_icpp(2)%alpha_rho(1)'   : awl2 * rho0wl2,           	
    'patch_icpp(2)%alpha_rho(2)'   : awv2 * rho0wv2,
    'patch_icpp(2)%alpha(1)'       : awl2,   	
    'patch_icpp(2)%alpha(2)'       : awv2,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    # ==========================================================
    # patch 3: with velocity ============================================
    'patch_icpp(3)%geometry'       : 3,     
    'patch_icpp(3)%x_centroid'     : (1.2E+00 + xe) * 2 / 4,     
    'patch_icpp(3)%y_centroid'     : ( ye + yb ) * 2 / 4,       
    'patch_icpp(3)%length_x'       : ( xe - 1.2E+00 ),       
    'patch_icpp(3)%length_y'       : ( ye - yb ) * 2 / 2,       
    'patch_icpp(3)%vel(1)'         : -200.0E00,   	
    'patch_icpp(3)%vel(2)'         : 0.0E00,   	
    'patch_icpp(3)%pres'           : p03,    	
    'patch_icpp(3)%alpha_rho(1)'   : awl3 * rho0wl3,           	
    'patch_icpp(3)%alpha_rho(2)'   : awv3 * rho0wv3,
    'patch_icpp(3)%alpha(1)'       : awl3,   	
    'patch_icpp(3)%alpha(2)'       : awv3,
    'patch_icpp(3)%alter_patch(1)' : 'T',       
    # Fluids Physical Parameters ===============================
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
    # ==========================================================
}))
