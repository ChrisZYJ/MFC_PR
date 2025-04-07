#:def Hardcoded2DVariables()

    real(wp) :: eps
    real(wp) :: r, rmax, gam, umax, p0
    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, intL, alph
    real(wp) :: alpha1, alpha2, alpha_rho1, alpha_rho2, y0, theta, f

    eps = 1e-9_wp

#:enddef

#:def Hardcoded2D()

    select case (patch_icpp(patch_id)%hcid) ! 2D_hardcoded_ic example case

    case (200)
        if (y_cc(j) <= (-x_cc(i)**3 + 1)**(1._wp/3._wp)) then
            ! Volume Fractions
            q_prim_vf(advxb)%sf(i, j, 0) = eps
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - eps
            ! Denssities
            q_prim_vf(contxb)%sf(i, j, 0) = eps*1000._wp
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - eps)*1._wp
            ! Pressure
            q_prim_vf(E_idx)%sf(i, j, 0) = 1000._wp
        end if
    case (202) ! Gresho vortex (Gouasmi et al 2022 JCP)
        r = ((x_cc(i) - 0.5_wp)**2 + (y_cc(j) - 0.5_wp)**2)**0.5_wp
        rmax = 0.2_wp

        gam = 1._wp + 1._wp/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1._wp/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5_wp)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5_wp)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5_wp)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2._wp)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2._wp + 4*(1 - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0._wp
            q_prim_vf(momxe)%sf(i, j, 0) = 0._wp
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*(-2 + 4*log(2._wp))
        end if
    case (203) ! Gresho vortex (Gouasmi et al 2022 JCP) with density correction
        r = ((x_cc(i) - 0.5_wp)**2._wp + (y_cc(j) - 0.5_wp)**2)**0.5_wp
        rmax = 0.2_wp

        gam = 1._wp + 1._wp/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1._wp/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5_wp)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5_wp)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5_wp)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2._wp/2._wp)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2._wp + 4._wp*(1._wp - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0._wp
            q_prim_vf(momxe)%sf(i, j, 0) = 0._wp
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2._wp*(-2._wp + 4*log(2._wp))
        end if

        q_prim_vf(contxb)%sf(i, j, 0) = q_prim_vf(E_idx)%sf(i, j, 0)**(1._wp/gam)
    case (204) ! Rayleigh-Taylor instability
        rhoH = 3._wp
        rhoL = 1._wp
        pRef = 1.e5_wp
        pInt = pRef
        h = 0.7_wp
        lam = 0.2_wp
        wl = 2._wp*pi/lam
        amp = 0.05_wp/wl

        intH = amp*sin(2._wp*pi*x_cc(i)/lam - pi/2._wp) + h

        alph = 0.5_wp*(1._wp + tanh((y_cc(j) - intH)/2.5e-3_wp))

        if (alph < eps) alph = eps
        if (alph > 1._wp - eps) alph = 1._wp - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, 0) = pref + rhoH*9.81_wp*(1.2_wp - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - alph)*rhoL
            pInt = pref + rhoH*9.81_wp*(1.2_wp - intH)
            q_prim_vf(E_idx)%sf(i, j, 0) = pInt + rhoL*9.81_wp*(intH - y_cc(j))
        end if

    case (205) ! 2D lung wave interaction problem
        h = 0.0           !non dim origin y
        lam = 1.0         !non dim lambda
        amp = patch_icpp(patch_id)%a(2)         !to be changed later!       !non dim amplitude

        intH = amp*sin(2*pi*x_cc(i)/lam - pi/2) + h

        if (y_cc(j) > intH) then
            q_prim_vf(contxb)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, 0) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, 0) = patch_icpp(1)%alpha(2)
        end if

    case (206) ! 2D lung wave interaction problem - horizontal domain
        h = 0.0           !non dim origin y
        lam = 1.0         !non dim lambda
        amp = patch_icpp(patch_id)%a(2)

        intL = amp*sin(2*pi*y_cc(j)/lam - pi/2) + h

        if (x_cc(i) > intL) then        !this is the liquid
            q_prim_vf(contxb)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, 0) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, 0) = patch_icpp(1)%alpha(2)
        end if

    case (281) ! Triangle Pointing Left
        alpha1 = 1.E-08
        alpha2 = 1-1.E-08
        alpha_rho1 = 1100*1.E-08
        alpha_rho2 = 1100*(1-1.E-08)

        ! fluid 1 if below the line y = 0.5*(x-0.04) and above the line y = -0.5*(x-0.04)
        if (y_cc(j) < 0.5*(x_cc(i)-0.04) .and. y_cc(j) > -0.5*(x_cc(i)-0.04)) then
            q_prim_vf(advxb)%sf(i, j, 0) = alpha1
            q_prim_vf(advxe)%sf(i, j, 0) = alpha2
            q_prim_vf(contxb)%sf(i, j, 0) = alpha_rho1
            q_prim_vf(contxe)%sf(i, j, 0) = alpha_rho2
        else
            q_prim_vf(advxb)%sf(i, j, 0) = alpha2
            q_prim_vf(advxe)%sf(i, j, 0) = alpha1
            q_prim_vf(contxb)%sf(i, j, 0) = alpha_rho2
            q_prim_vf(contxe)%sf(i, j, 0) = alpha_rho1
        end if

    case (282) ! Triangle Pointing Right
        alpha1 = 1.E-08
        alpha2 = 1-1.E-08
        alpha_rho1 = 1100*1.E-08
        alpha_rho2 = 1100*(1-1.E-08)

        ! fluid 1 if below the line y = -0.5*(x-0.06) and above the line y = 0.5*(x-0.06)
        if (y_cc(j) < -0.5*(x_cc(i)-0.06) .and. y_cc(j) > 0.5*(x_cc(i)-0.06)) then
            q_prim_vf(advxb)%sf(i, j, 0) = alpha1
            q_prim_vf(advxe)%sf(i, j, 0) = alpha2
            q_prim_vf(contxb)%sf(i, j, 0) = alpha_rho1
            q_prim_vf(contxe)%sf(i, j, 0) = alpha_rho2
        else
            q_prim_vf(advxb)%sf(i, j, 0) = alpha2
            q_prim_vf(advxe)%sf(i, j, 0) = alpha1
            q_prim_vf(contxb)%sf(i, j, 0) = alpha_rho2
            q_prim_vf(contxe)%sf(i, j, 0) = alpha_rho1
        end if

    case (283) ! Parabola (right) with focal point at (0.05, 0)
        alpha1 = 1.E-08
        alpha2 = 1-1.E-08
        alpha_rho1 = 1100*1.E-08
        alpha_rho2 = 1100*(1-1.E-08)

        ! fluid 1 if right of the parabola x = 0.05 + (y0**2 - y_cc(j)**2)/(2.0*y0)
        ! where +-y0 are the y-intercepts
        y0 = 0.02
        if (x_cc(i) > 0.05 + (y0**2 - y_cc(j)**2)/(2.0*y0)) then
            q_prim_vf(advxb)%sf(i, j, 0) = alpha1
            q_prim_vf(advxe)%sf(i, j, 0) = alpha2
            q_prim_vf(contxb)%sf(i, j, 0) = alpha_rho1
            q_prim_vf(contxe)%sf(i, j, 0) = alpha_rho2
        else
            q_prim_vf(advxb)%sf(i, j, 0) = alpha2
            q_prim_vf(advxe)%sf(i, j, 0) = alpha1
            q_prim_vf(contxb)%sf(i, j, 0) = alpha_rho2
            q_prim_vf(contxe)%sf(i, j, 0) = alpha_rho1
        end if

    case (284) ! Bouba-like Shape
        ! Define fluid properties for two fluids
        alpha1 = 1.E-08
        alpha2 = 1-1.E-08
        alpha_rho1 = 1100*1.E-08
        alpha_rho2 = 1100*(1-1.E-08)

        ! Calculate polar angle theta using arctan2 for correct quadrant
        theta = atan2(y_cc(j), x_cc(i)-0.05)

        ! Compute the implicit function F(x, y) for the bouba shape
        F = sqrt((x_cc(i)-0.05)**2 + y_cc(j)**2) / 0.01 - &
            (1.0 + 0.3 * cos(5.0 * theta + 0.0) + &
             0.15 * cos(7.0 * theta + pi/3) + &
             0.25 * cos(11.0 * theta + pi/2))

        ! Determine if the point is inside or outside the bouba shape
        if (F < 0.0) then
            ! Inside the bouba shape: Assign Fluid 1 properties
            q_prim_vf(advxb)%sf(i, j, 0) = alpha1
            q_prim_vf(advxe)%sf(i, j, 0) = alpha2
            q_prim_vf(contxb)%sf(i, j, 0) = alpha_rho1
            q_prim_vf(contxe)%sf(i, j, 0) = alpha_rho2
        else
            ! Outside the bouba shape: Assign Fluid 2 properties
            q_prim_vf(advxb)%sf(i, j, 0) = alpha2
            q_prim_vf(advxe)%sf(i, j, 0) = alpha1
            q_prim_vf(contxb)%sf(i, j, 0) = alpha_rho2
            q_prim_vf(contxe)%sf(i, j, 0) = alpha_rho1
        end if

    case (285) ! Bouba-like Shape (more round and off-centered)
        ! Define fluid properties for two fluids
        alpha1 = 1.E-08
        alpha2 = 1-1.E-08
        alpha_rho1 = 1100*1.E-08
        alpha_rho2 = 1100*(1-1.E-08)

        ! Calculate polar angle theta using arctan2 for correct quadrant
        theta = atan2((y_cc(j)), x_cc(i)-0.05)

        ! Compute the implicit function F(x, y) for the bouba shape
        F = sqrt((x_cc(i)-0.05)**2 + y_cc(j)**2) / 0.01 - &
            (1.0 + 0.1 * cos(5.0 * theta + 0.0) + &
             0.02 * cos(7.0 * theta + pi/3) + &
             0.05 * cos(11.0 * theta + pi/2) + &
             0.8 * cos(theta - pi/10))

        ! Determine if the point is inside or outside the bouba shape
        if (F < 0.0) then
            ! Inside the bouba shape: Assign Fluid 1 properties
            q_prim_vf(advxb)%sf(i, j, 0) = alpha1
            q_prim_vf(advxe)%sf(i, j, 0) = alpha2
            q_prim_vf(contxb)%sf(i, j, 0) = alpha_rho1
            q_prim_vf(contxe)%sf(i, j, 0) = alpha_rho2
        else
            ! Outside the bouba shape: Assign Fluid 2 properties
            q_prim_vf(advxb)%sf(i, j, 0) = alpha2
            q_prim_vf(advxe)%sf(i, j, 0) = alpha1
            q_prim_vf(contxb)%sf(i, j, 0) = alpha_rho2
            q_prim_vf(contxe)%sf(i, j, 0) = alpha_rho1
        end if

    case default
        if (proc_rank == 0) then
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch "//trim(iStr))
        end if

    end select

#:enddef
