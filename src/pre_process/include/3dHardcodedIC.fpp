#:def Hardcoded3DVariables()
    ! Place any declaration of intermediate variables here

    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, alph

    real(wp) :: eps
    real(wp) :: r, x_center, y_center, z_center, radius
    real(wp) :: delta_x, delta_y, delta_z, theta, phi, r_irreg

    eps = 1e-9_wp
#:enddef

#:def Hardcoded3D()

    select case (patch_icpp(patch_id)%hcid)
    case (300) ! Rayleigh-Taylor instability
        rhoH = 3._wp
        rhoL = 1._wp
        pRef = 1.e5_wp
        pInt = pRef
        h = 0.7_wp
        lam = 0.2_wp
        wl = 2._wp*pi/lam
        amp = 0.025_wp/wl

        intH = amp*(sin(2._wp*pi*x_cc(i)/lam - pi/2._wp) + sin(2._wp*pi*z_cc(k)/lam - pi/2._wp)) + h

        alph = 5e-1_wp*(1._wp + tanh((y_cc(j) - intH)/2.5e-3_wp))

        if (alph < eps) alph = eps
        if (alph > 1._wp - eps) alph = 1._wp - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1._wp - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, k) = pref + rhoH*9.81_wp*(1.2_wp - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, k) = alph
            q_prim_vf(advxe)%sf(i, j, k) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, k) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, k) = (1._wp - alph)*rhoL
            pInt = pref + rhoH*9.81_wp*(1.2_wp - intH)
            q_prim_vf(E_idx)%sf(i, j, k) = pInt + rhoL*9.81_wp*(intH - y_cc(j))
        end if

    case (301) ! (3D lung geometry in X direction, |sin(*)+sin(*)|)
        h = 0.0_wp
        lam = 1.0_wp
        amp = patch_icpp(patch_id)%a(2)
        intH = amp*abs((sin(2*pi*y_cc(j)/lam - pi/2) + sin(2*pi*z_cc(k)/lam - pi/2)) + h)
        if (x_cc(i) > intH) then
            q_prim_vf(contxb)%sf(i, j, k) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, k) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, k) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, k) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, k) = patch_icpp(1)%alpha(2)
        end if

    case (380)
        ! INFO:
        !   x_centroid: 0.05
        !   y_centroid: 0.0
        !   z_centroid: 0.0
        !   radius    : 0.03
        !
        !   Within the circle:
        !     alpha_rho(1) : 1100*1.E-08
        !     alpha_rho(2) : 1100*(1 - 1.E-08)
        !     alpha(1)     : 1.E-08
        !     alpha(2)     : 1.0 - 1.E-08
        !
        !   Outside the circle:
        !     alpha_rho(1) : 1100*(1 - 1.E-08)
        !     alpha_rho(2) : 1100*1.E-08
        !     alpha(1)     : 1.0 - 1.E-08
        !     alpha(2)     : 1.E-08

        x_center = 0.05_wp
        y_center = 0.0_wp
        z_center = 0.0_wp
        radius   = 0.02_wp

        r = sqrt((x_cc(i) - x_center)**2 + (y_cc(j) - y_center)**2 + (z_cc(k) - z_center)**2)

        if (r < radius) then
            ! Within the circle:
            q_prim_vf(contxb)%sf(i, j, k) = 1100*1.E-08_wp
            q_prim_vf(contxe)%sf(i, j, k) = 1100*(1.0_wp - 1.E-08_wp)
            q_prim_vf(advxb)%sf(i, j, k) = 1.E-08_wp
            q_prim_vf(advxe)%sf(i, j, k) = 1.0_wp - 1.E-08_wp
        else
            ! Outside the circle:
            q_prim_vf(contxb)%sf(i, j, k) = 1100*(1.0_wp - 1.E-08_wp)
            q_prim_vf(contxe)%sf(i, j, k) = 1100*1.E-08_wp
            q_prim_vf(advxb)%sf(i, j, k) = 1.0_wp - 1.E-08_wp
            q_prim_vf(advxe)%sf(i, j, k) = 1.E-08_wp
        end if

    case (381)
        ! INFO:
        !   x_centroid: 0.05
        !   y_centroid: 0.0
        !   z_centroid: 0.0
        !   The irregular solid is defined via spherical harmonic-like modulations.
        !   Maximum radius is capped at 0.035, with a baseline of 0.03 and deformations added.
        !
        !   Within the shape:
        !     alpha_rho(1) : 1100*1.E-08
        !     alpha_rho(2) : 1100*(1 - 1.E-08)
        !     alpha(1)     : 1.E-08
        !     alpha(2)     : 1.0 - 1.E-08
        !
        !   Outside the shape:
        !     alpha_rho(1) : 1100*(1 - 1.E-08)
        !     alpha_rho(2) : 1100*1.E-08
        !     alpha(1)     : 1.0 - 1.E-08
        !     alpha(2)     : 1.E-08

        ! Set the center of the irregular shape:
        x_center = 0.05_wp
        y_center = 0.0_wp
        z_center = 0.0_wp

        ! Compute the offset from the center:
        delta_x = x_cc(i) - x_center
        delta_y = y_cc(j) - y_center
        delta_z = z_cc(k) - z_center
        r = sqrt(delta_x**2 + delta_y**2 + delta_z**2)

        ! Compute spherical angles relative to the center
        if (r > 1e-12_wp) then
            theta = acos(delta_z/r)
            phi   = atan2(delta_y, delta_x)
        else
            theta = 0.0_wp
            phi   = 0.0_wp
        endif

        ! Define an irregular radius function using a combination of harmonic modulations.
        ! Baseline radius is 0.03_wp; modulation coefficients yield deformations and spikes.
        r_irreg = 0.02_wp * ( 1.0_wp &
                 + 0.07_wp*sin(6.7_wp*theta + 3.2_wp)*sin(6.5_wp*phi + 0.16_wp) &
                 - 0.09_wp*cos(4.1_wp*theta + 5.8_wp)*cos(3.1_wp*phi + 0.33_wp) &
                 + 0.05_wp*sin(2.2_wp*theta + 1.1_wp)*sin(5.0_wp*phi) &
                 - 0.13_wp*sin(6.1_wp*theta)*cos(2.2_wp*phi) &
                 + 0.19_wp*sin(4.2_wp*theta + 1.1_wp*phi + 9.0) &
                 + 0.05_wp*cos(theta + 8.9_wp) )
        ! Ensure the radius does not exceed 0.035_wp:
        if (r_irreg > 0.03_wp) then
            r_irreg = 0.03_wp
        endif

        ! Determine if the current grid point is inside the irregular solid:
        if (r < r_irreg) then
            ! Inside the irregular solid:
            q_prim_vf(contxb)%sf(i, j, k) = 1100*1.E-08_wp
            q_prim_vf(contxe)%sf(i, j, k) = 1100*(1.0_wp - 1.E-08_wp)
            q_prim_vf(advxb)%sf(i, j, k) = 1.E-08_wp
            q_prim_vf(advxe)%sf(i, j, k) = 1.0_wp - 1.E-08_wp
        else
            ! Outside the irregular solid:
            q_prim_vf(contxb)%sf(i, j, k) = 1100*(1.0_wp - 1.E-08_wp)
            q_prim_vf(contxe)%sf(i, j, k) = 1100*1.E-08_wp
            q_prim_vf(advxb)%sf(i, j, k) = 1.0_wp - 1.E-08_wp
            q_prim_vf(advxe)%sf(i, j, k) = 1.E-08_wp
        endif

        ! Put your variable assignments here
    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
    end select

#:enddef
