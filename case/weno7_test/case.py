#!/usr/bin/env python3

import json

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : 0.E+00,
    'x_domain%end'                 : 1.E+00,
    'm'                            : 36,
    'n'                            : 0,
    'p'                            : 0,
    'dt'                           : 0.003,
    't_step_start'                 : 0,
    't_step_stop'                  : 50,
    't_step_save'                  : 1,
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 2,
    'model_eqns'                   : 2,
    'alt_soundspeed'               : 'F',
    'num_fluids'                   : 1,
    'mpp_lim'                      : 'F',
    'mixture_err'                  : 'F',
    'time_stepper'                 : 1,
    'weno_order'                   : 7,
    'weno_eps'                     : 1.E-40,
    'wenoz'                        : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -3,
    'bc_x%end'                     : -3,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 2,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'rho_wrt'                      :'T',
    'parallel_io'                  :'F',
    # ==========================================================================
                                                            
    # Patch 1 L ================================================================
    'patch_icpp(1)%geometry'       : 1,
    'patch_icpp(1)%x_centroid'     : 0.25,
    'patch_icpp(1)%length_x'       : 0.5,
    'patch_icpp(1)%vel(1)'         : 0.698,
    'patch_icpp(1)%pres'           : 3.528,
    'patch_icpp(1)%alpha_rho(1)'   : 0.445E+00,
    'patch_icpp(1)%alpha(1)'       : 1.,
    # ==========================================================================

    # Patch 2 R ================================================================
    'patch_icpp(2)%geometry'       : 1,
    'patch_icpp(2)%x_centroid'     : 0.75,
    'patch_icpp(2)%length_x'       : 0.5,
    'patch_icpp(2)%vel(1)'         : 0.0,
    'patch_icpp(2)%pres'           : 0.571,
    'patch_icpp(2)%alpha_rho(1)'   : 0.5E+00,
    'patch_icpp(2)%alpha(1)'       : 1.,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.E+00/(1.4-1.E+00),
    'fluid_pp(1)%pi_inf'           : 0.0,
    # ==========================================================================
}))

# ==============================================================================