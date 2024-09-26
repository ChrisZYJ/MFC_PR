#!/usr/bin/env python3

import math
import json

# Define non-dimensional parameters based on Table 2
Ma = 0.1       # Mach number
Re = 10000     # Reynolds number
Pr = 0.73      # Prandtl number
mygamma = 1.4  # Specific heat ratio (renamed to avoid duplication)

# Dynamic viscosity based on Reynolds number
Mu = 1.0 / Re  # Mu = 1/Re = 1/10000 = 0.0001

# Pressure based on Mach number
p = 1.0 / (mygamma * Ma**2)  # p = 1 / (1.4 * 0.01) = 71.42857142857143

# Grid resolution for theta=120
m = 160
n = 160

# Configuring case dictionary
case_dict = {
    # Logistics ================================================
    'run_time_info'                : 'T',
    # ==========================================================

    # Computational Domain Parameters ==========================
    # Unit square domain for Double Periodic Shear Layer
    'x_domain%beg'                 : 0.0, 
    'x_domain%end'                 : 1.0,
    'y_domain%beg'                 : 0.0,
    'y_domain%end'                 : 1.0,
    'm'                            : m,  # Grid resolution
    'n'                            : n,
    'p'                            : 0,
    'dt'                           : 2.E-04,
    't_step_start'                 : 0,
    't_step_stop'                  : 5000,
    't_step_save'                  : 10,
    # ==========================================================

    # Simulation Algorithm Parameters ==========================
    'num_patches'                  : 1,
    'model_eqns'                   : 2,  # Assuming Euler equations; adjust if different
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 1,
    'adv_alphan'                   : 'T',
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'T',
    'time_stepper'                 : 3,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    # 'mapped_weno'                  : 'T',
    'teno'                         : 'T',
    'teno_CT'                      : 1.E-6,
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'weno_Re_flux'                 : 'T',
    'weno_avg'                     : 'T',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -1,  # Periodic boundaries
    'bc_x%end'                     : -1,
    'bc_y%beg'                     : -1,
    'bc_y%end'                     : -1,
    # ==========================================================

    # Formatted Database Files Structure Parameters ============
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                : 'T',
    'parallel_io'                  : 'T',
    'fd_order'                     : 2,
    'omega_wrt(3)'                 : 'T',
    # ==========================================================

    # Patch 1: Double Periodic Shear Layer ======================
    'patch_icpp(1)%hcid'           : 205,                           # Set to Case 205
    'patch_icpp(1)%geometry'       : 7,                             # 2D hard-coded
    'patch_icpp(1)%x_centroid'     : 0.5,                           # Centroid at (0.5, 0.5) for periodic domain
    'patch_icpp(1)%y_centroid'     : 0.5,
    'patch_icpp(1)%length_x'       : 1.0,                           # Domain length in x
    'patch_icpp(1)%length_y'       : 1.0,                           # Domain length in y
    'patch_icpp(1)%vel(1)'         : 0.0,   # Dummy
    'patch_icpp(1)%vel(2)'         : 0.0,   # Dummy
    'patch_icpp(1)%pres'           : p,                             # Pressure based on Mach number
    # Single-phase flow parameters
    'patch_icpp(1)%alpha_rho(1)'   : 1.0,                           # rhoH = rhoL = rho = 1.0
    'patch_icpp(1)%alpha(1)'       : 1.0,                           # Volume fraction alpha = 1.0
    # ==========================================================

    # Fluids Physical Parameters ===============================
    'fluid_pp(1)%gamma'            : 1.E+00/(mygamma-1.E+00),
    'fluid_pp(1)%pi_inf'           : 0.0,
    'fluid_pp(1)%Re(1)'            : Re,                            # Shear viscosity: Re
    # ==========================================================
}

# Output the case dictionary as JSON
print(json.dumps(case_dict, indent=4))
