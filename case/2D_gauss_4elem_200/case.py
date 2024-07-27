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
                    'y_domain%beg'                 : -50.E-03,                 \
                    'y_domain%end'                 : 50.E-03,                  \
                    'm'                            : 199,                      \
                    'n'                            : 199,                      \
                    'p'                            : 0,                        \
                    'dt'                           : 2E-7,                     \
                    't_step_start'                 : 0,                        \
                    't_step_stop'                  : 600,                      \
                    't_step_save'                  : 4,                        \
                    # ==========================================================
                                                                               \
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 2,                        \
                    'model_eqns'                   : 2,                        \
                    'alt_soundspeed'               : 'F',                      \
                    'num_fluids'                   : 2,                        \
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
                    'bc_y%beg'                     : -8,                       \
                    'bc_y%end'                     : -8,                       \
                    'weno_Re_flux'                 : 'F',                      \
                    # ==========================================================

                    # Turning on Hypoelasticity ================================
                    'hypoelasticity'               : 'T',                      \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    # 'format'                       : 1,                        \
                    # 'prim_vars_wrt'                :'T',                       \

                    'format'                       : 2,                        \
                    'output_partial_domain'        :'T',                       \
                    'x_output%beg'                 : 0.047E+00,                \
                    'x_output%end'                 : 0.073E+00,                \
                    'y_output%beg'                 :-0.003E+00,                \
                    'y_output%end'                 : 0.003E+00,                \

                    'precision'                    : 2,                        \
                    'vel_wrt(1)'                   :'T',                       \
                    'vel_wrt(2)'                   :'T',                       \
                    'pres_wrt'                     :'T',                       \
                    'parallel_io'                  :'T',                       \
                    'fd_order'                     : 1,                        \
                    # ==========================================================
                                                                                        
                    # Patch 1 Liquid ===========================================
                    'patch_icpp(1)%geometry'       : 3,                         \
                    'patch_icpp(1)%x_centroid'     : 0.05,                      \
                    'patch_icpp(1)%y_centroid'     : 0.0,                       \
                    'patch_icpp(1)%length_x'       : 0.1,                       \
                    'patch_icpp(1)%length_y'       : 0.1,                       \
                    'patch_icpp(1)%vel(1)'         : 0.0,                       \
                    'patch_icpp(1)%vel(2)'         : 0.0,                       \
                    'patch_icpp(1)%pres'           : 1E+05,                     \
                    'patch_icpp(1)%alpha_rho(1)'   : 1100*(1-1.E-08),           \
                    'patch_icpp(1)%alpha_rho(2)'   : 1100*1.E-08,               \
                    'patch_icpp(1)%alpha(1)'       : 1.0-1.E-08,                \
                    'patch_icpp(1)%alpha(2)'       : 1.E-08,                    \
                    'patch_icpp(1)%tau_e(1)'       : 0.0,                       \
                    'patch_icpp(1)%tau_e(2)'       : 0.0,                       \
                    'patch_icpp(1)%tau_e(3)'       : 0.0,                       \
                    # ==========================================================

                    # Patch 2 Stone ============================================
                    'patch_icpp(2)%geometry'       : 3,                         \
                    'patch_icpp(2)%x_centroid'     : 0.06,                      \
                    'patch_icpp(2)%y_centroid'     : 0.0,                       \
                    'patch_icpp(2)%length_x'       : 0.026,                     \
                    'patch_icpp(2)%length_y'       : 0.006,                     \
                    'patch_icpp(2)%alter_patch(1)' : 'T',                       \
                    'patch_icpp(2)%vel(1)'         : 0,                         \
                    'patch_icpp(2)%vel(2)'         : 0,                         \
                    'patch_icpp(2)%pres'           : 1E+05,                     \
                    'patch_icpp(2)%alpha_rho(1)'   : 1100*1.E-08,               \
                    'patch_icpp(2)%alpha_rho(2)'   : 1100*(1-1.E-08),           \
                    'patch_icpp(2)%alpha(1)'       : 1.E-08,                    \
                    'patch_icpp(2)%alpha(2)'       : 1-1.E-08,                  \
                    'patch_icpp(2)%tau_e(1)'       : 0.0,                       \
                    'patch_icpp(2)%tau_e(2)'       : 0.0,                       \
                    'patch_icpp(2)%tau_e(3)'       : 0.0,                       \
                    # ==========================================================

                    # ==========================================================
                    'acoustic_source'                   : 'T',
                    'num_source'                        : 1,
                    'acoustic(1)%support'               : 9,
                    'acoustic(1)%loc(1)'                : 0.006,
                    'acoustic(1)%loc(2)'                : 0,
                    'acoustic(1)%pulse'                 : 2,
                    'acoustic(1)%npulse'                : 1,
                    'acoustic(1)%mag'                   : 10,
                    'acoustic(1)%gauss_sigma_time'      : 4E-7,
                    'acoustic(1)%foc_length'            : 0.054,
                    'acoustic(1)%aperture'              : 0.08,
                    'acoustic(1)%num_elements'          : 4,
                    'acoustic(1)%element_on'            : 0,
                    'acoustic(1)%element_spacing_angle' : math.pi/180*5,
                    'acoustic(1)%delay'                 : 16E-7,
                    # ==========================================================

                    # Fluids Physical Parameters ===============================
                    'fluid_pp(1)%gamma'            : 1.E+00/(4.4E+00-1.E+00),               \
                    'fluid_pp(1)%pi_inf'           : 4.4E+00*5.57E+08/(4.4E+00 - 1.E+00),   \
                    'fluid_pp(1)%G'                : 0.,                                    \
                    'fluid_pp(2)%gamma'            : 1.E+00/(4.4E+00-1.E+00),               \
                    'fluid_pp(2)%pi_inf'           : 4.4E+00*9.2940E+08/(4.4E+00 - 1.E+00), \
                    'fluid_pp(2)%G'                : 1.8447E+09,                            \
                    # ==========================================================
    }))

