#!/usr/bin/python3

# This input file is used for 2D impulse response with a Gaussian of width 1e-7.

import json, math

# Configuring case dictionary
print(json.dumps({
                    # Logistics ================================================
                    # 'case_dir'                     : '\'.\'',                  \
                    'run_time_info'                : 'T',                      \
                    # ==========================================================
                                                                               \
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 : 0.E+00,                   \
                    'x_domain%end'                 : 100.E-03,                 \
                    'm'                            : 49,                       \
                    'n'                            : 0,                        \
                    'p'                            : 0,                        \
                    'dt'                           : 8E-7,                     \
                    't_step_start'                 : 0,                        \
                    't_step_stop'                  : 150,                      \
                    't_step_save'                  : 1,                        \
                    # ==========================================================
                                                                               \
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 1,                        \
                    'model_eqns'                   : 2,                        \
                    'alt_soundspeed'               : 'F',                      \
                    'num_fluids'                   : 1,                        \
                    'adv_alphan'                   : 'T',                      \
                    'mpp_lim'                      : 'F',                      \
                    'mixture_err'                  : 'F',                      \
                    'time_stepper'                 : 3,                        \
                    # 'weno_vars'                    : 2,                        \
                    'weno_order'                   : 5,                        \
                    'weno_eps'                     : 1.E-40,                   \
                    'mapped_weno'                  : 'F',                      \
                    'teno'                         : 'T',                      \
                    'teno_CT'                      : 1.E-6,                    \
                    'null_weights'                 : 'F',                      \
                    'mp_weno'                      : 'F',                      \
                    'riemann_solver'               : 1,                        \
                    'wave_speeds'                  : 1,                        \
                    'avg_state'                    : 2,                        \
                    'bc_x%beg'                     : -8,                       \
                    'bc_x%end'                     : -8,                       \
                    'weno_Re_flux'                 : 'F',                      \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 2,                        \
                    'precision'                    : 2,                        \
                    'prim_vars_wrt'                :'T',                       \
                    'parallel_io'                  :'T',                       \
                    'fd_order'                     : 1,                        \
                    # ==========================================================
                                                                                        
                    # Patch 1 Liquid ===========================================
                    'patch_icpp(1)%geometry'       : 1,                         \
                    'patch_icpp(1)%x_centroid'     : 0.05,                      \
                    'patch_icpp(1)%length_x'       : 0.1,                       \
                    'patch_icpp(1)%vel(1)'         : 0.0,                       \
                    'patch_icpp(1)%pres'           : 1E+05,                     \
                    'patch_icpp(1)%alpha_rho(1)'   : 1100,                      \
                    'patch_icpp(1)%alpha(1)'       : 1.0,                       \
                    # ==========================================================

                    # ==========================================================
                    'acoustic_source'                   : 'T',
                    'num_source'                        : 1,
                    'acoustic(1)%support'               : 1,
                    'acoustic(1)%loc(1)'                : 0.01,
                    'acoustic(1)%pulse'                 : 2,
                    'acoustic(1)%npulse'                : 1,
                    'acoustic(1)%mag'                   : 10,
                    'acoustic(1)%gauss_sigma_time'      : 4E-7,
                    'acoustic(1)%dir'                   : 1.,
                    'acoustic(1)%delay'                 : 16E-7,
                    # ==========================================================

                    # Fluids Physical Parameters ===============================
                    'fluid_pp(1)%gamma'            : 1.E+00/(4.4E+00-1.E+00),               \
                    'fluid_pp(1)%pi_inf'           : 4.4E+00*5.57E+08/(4.4E+00 - 1.E+00),   \
                    'fluid_pp(1)%G'                : 0.,                                    \
                    # ==========================================================
    }))

