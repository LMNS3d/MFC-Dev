#!/usr/bin/env python3
import math, json
## 1 FOR BACKGROUND, 2 FOR SHOKED STATE, 3 FOR WATER REGION (WHEN NEEDED)
# Pressure
p01 = 1.E+05
p02 = 1.E+05
p03 = 1.E+05

## FLUID PROPERTIES FOR EACH PATCH

## liquid dodecane

# pi infty
piwl = 4.0E+08
# qv
qvwl = -775269
# qv'
qvpwl = 0.0E0
# cv
cvwl = 1077.7
# gamma
gamwl = 2.35
# cp
cpwl = cvwl*gamwl

# density
rho0wl1 = 458.338E+00
rho0wl2 = 458.338E+00
rho0wl3 = 458.338E+00

# volume fraction
awl1 = 1.0E+00 - 1.0E-06
awl2 = 1.0E-06
awl3 = 1.0E+00 - 1.0E-06

## Vapor dodecane

# pi infinity
piwv = 0
# qv
qvwv = -237547
# qv'
qvpwv = -24400
# cv
cvwv = 1956.45
# gamma
gamwv = 1.025
# cp
cpwv = cvwv*gamwv

# density
rho0wv1 = 3.408E+0
rho0wv2 = 3.408E+0
rho0wv3 = 3.408E+0

# volume fraction
awv1 = 1.0E-06
awv2 = 1.0E+00 - 1.0E-06
awv3 = 1.0E-06

# discretization parameters for the x direction
# number of elements
Nx      = 199

# domain boundaries
xb = 0.00
xe = 1.00

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
    'dt'           : 1.E-06,        
    't_step_start' : 0,         
    't_step_stop'  : 1500,        
    't_step_save'  : 50,         
    # ==========================================================
    # Simulation Algorithm Parameters ==========================
    'num_patches'  : 2,        
    'model_eqns'   : 3,        
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
    'bc_x%beg'     : -16,    
    'bc_x%vb1'     : 100,
    'bc_x%end'     : -2,       
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
    'patch_icpp(1)%vel(1)'         : 0.0E+00,   
    'patch_icpp(1)%vel(2)'         : 0.0E+00,   
    'patch_icpp(1)%pres'           : p01,  	
    'patch_icpp(1)%alpha_rho(1)'   : awl1 * rho0wl1,           	
    'patch_icpp(1)%alpha_rho(2)'   : awv1 * rho0wv1,
    'patch_icpp(1)%alpha(1)'       : awl1,   	
    'patch_icpp(1)%alpha(2)'       : awv1,
    # ==========================================================
    # Patch 2 - vapor bubble ==================================
    'patch_icpp(2)%geometry'       : 2,     	
    'patch_icpp(2)%x_centroid'     : ( xe + xb ) * 2 / 4,       
    'patch_icpp(2)%y_centroid'     : ( ye + yb ) * 2 / 4,       
    'patch_icpp(2)%radius'         : 0.2E00,
    'patch_icpp(2)%vel(1)'         : 0.0E00,   	
    'patch_icpp(2)%vel(2)'         : 0.0E00,   	    
    'patch_icpp(2)%pres'           : p02,    	
    'patch_icpp(2)%alpha_rho(1)'   : awl2 * rho0wl2,           	
    'patch_icpp(2)%alpha_rho(2)'   : awv2 * rho0wv2,
    'patch_icpp(2)%alpha(1)'       : awl2,   	
    'patch_icpp(2)%alpha(2)'       : awv2,
    'patch_icpp(2)%alter_patch(1)' : 'T',
    # ==========================================================
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
