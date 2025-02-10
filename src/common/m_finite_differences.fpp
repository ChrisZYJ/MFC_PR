module m_finite_differences

    use m_global_parameters

    implicit none

contains

    subroutine s_compute_fd_divergence(div, fields, ix_s, iy_s, iz_s)

        type(scalar_field), intent(INOUT) :: div
        type(scalar_field), intent(IN) :: fields(1:3)
        type(int_bounds_info), intent(IN) :: ix_s, iy_s, iz_s

        integer :: x, y, z !< Generic loop iterators

        real(wp) :: divergence

        !$acc parallel loop collapse(3) private(divergence)
        do x = ix_s%beg, ix_s%end
            do y = iy_s%beg, iy_s%end
                do z = iz_s%beg, iz_s%end

                    if (x == ix_s%beg) then
                        divergence = (-3._wp*fields(1)%sf(x, y, z) + 4._wp*fields(1)%sf(x + 1, y, z) - fields(1)%sf(x + 2, y, z))/(x_cc(x + 2) - x_cc(x))
                    else if (x == ix_s%end) then
                        divergence = (+3._wp*fields(1)%sf(x, y, z) - 4._wp*fields(1)%sf(x - 1, y, z) + fields(1)%sf(x - 2, y, z))/(x_cc(x) - x_cc(x - 2))
                    else
                        divergence = (fields(1)%sf(x + 1, y, z) - fields(1)%sf(x - 1, y, z))/(x_cc(x + 1) - x_cc(x - 1))
                    end if

                    if (n > 0) then
                        if (y == iy_s%beg) then
                            divergence = divergence + (-3._wp*fields(2)%sf(x, y, z) + 4._wp*fields(2)%sf(x, y + 1, z) - fields(2)%sf(x, y + 2, z))/(y_cc(y + 2) - y_cc(y))
                        else if (y == iy_s%end) then
                            divergence = divergence + (+3._wp*fields(2)%sf(x, y, z) - 4._wp*fields(2)%sf(x, y - 1, z) + fields(2)%sf(x, y - 2, z))/(y_cc(y) - y_cc(y - 2))
                        else
                            divergence = divergence + (fields(2)%sf(x, y + 1, z) - fields(2)%sf(x, y - 1, z))/(y_cc(y + 1) - y_cc(y - 1))
                        end if
                    end if

                    if (p > 0) then
                        if (z == iz_s%beg) then
                            divergence = divergence + (-3._wp*fields(3)%sf(x, y, z) + 4._wp*fields(3)%sf(x, y, z + 1) - fields(3)%sf(x, y, z + 2))/(z_cc(z + 2) - z_cc(z))
                        else if (z == iz_s%end) then
                            divergence = divergence + (+3._wp*fields(3)%sf(x, y, z) - 4._wp*fields(3)%sf(x, y, z - 1) + fields(2)%sf(x, y, z - 2))/(z_cc(z) - z_cc(z - 2))
                        else
                            divergence = divergence + (fields(3)%sf(x, y, z + 1) - fields(3)%sf(x, y, z - 1))/(z_cc(z + 1) - z_cc(z - 1))
                        end if
                    end if

                    div%sf(x, y, z) = div%sf(x, y, z) + divergence

                end do
            end do
        end do

    end subroutine s_compute_fd_divergence

    !>  The purpose of this subroutine is to compute the finite-
    !!      difference coefficients for the centered schemes utilized
    !!      in computations of first order spatial derivatives in the
    !!      s-coordinate direction. The s-coordinate direction refers
    !!      to the x-, y- or z-coordinate direction, depending on the
    !!      subroutine's inputs. Note that coefficients of up to 4th
    !!      order accuracy are available.
    !!  @param q Number of cells in the s-coordinate direction
    !!  @param s_cc Locations of the cell-centers in the s-coordinate direction
    !!  @param fd_coeff_s Finite-diff. coefficients in the s-coordinate direction
    subroutine s_compute_finite_difference_coefficients(q, s_cc, fd_coeff_s, buff_size, &
                                                        fd_number_in, fd_order_in, offset_s)

        integer :: lB, lE !< loop bounds
        integer, intent(IN) :: q
        integer, intent(IN) :: buff_size, fd_number_in, fd_order_in
        type(int_bounds_info), optional, intent(IN) :: offset_s
        real(wp), allocatable, dimension(:, :), intent(INOUT) :: fd_coeff_s

        real(wp), &
            dimension(-buff_size:q + buff_size), &
            intent(IN) :: s_cc

        integer :: i !< Generic loop iterator
        real(wp) :: dx0

        if (present(offset_s)) then
            lB = -offset_s%beg
            lE = q + offset_s%end
        else
            lB = 0
            lE = q
        end if

#ifdef MFC_POST_PROCESS
        if (allocated(fd_coeff_s)) deallocate (fd_coeff_s)
        allocate (fd_coeff_s(-fd_number_in:fd_number_in, lb:lE))
#endif

        ! Computing the 1st order finite-difference coefficients
        if (fd_order_in == 1) then
            do i = lB, lE
                fd_coeff_s(-1, i) = 0._wp
                fd_coeff_s(0, i) = -1._wp/(s_cc(i + 1) - s_cc(i))
                fd_coeff_s(1, i) = -fd_coeff_s(0, i)
            end do

            ! Computing the 2nd order finite-difference coefficients
        elseif (fd_order_in == 2) then
            do i = lB, lE
                fd_coeff_s(-1, i) = -1._wp/(s_cc(i + 1) - s_cc(i - 1))
                fd_coeff_s(0, i) = 0._wp
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
            end do

            ! Computing the 4th order finite-difference coefficients
        else
            do i = lB, lE
                fd_coeff_s(-2, i) = 1._wp/(s_cc(i - 2) - 8._wp*s_cc(i - 1) - s_cc(i + 2) + 8._wp*s_cc(i + 1))
                fd_coeff_s(-1, i) = -8._wp*fd_coeff_s(-2, i)
                fd_coeff_s(0, i) = 0._wp
                fd_coeff_s(1, i) = -fd_coeff_s(-1, i)
                fd_coeff_s(2, i) = -fd_coeff_s(-2, i)
            end do

        end if

        print *, (s_cc(2)-s_cc(1))/(s_cc(1)-s_cc(0))
        if (abs((s_cc(2)-s_cc(1))/(s_cc(1)-s_cc(0)) - 1._wp) > 0.1_wp) then ! cyl_coord y-dir
            print *, 'DEBUG'
            ! ! For the uniform grid case (already given):
            ! fd_coeff_s(-1, i) = -1._wp/(s_cc(i+1)-s_cc(i-1))
            ! fd_coeff_s(0,  i) =  0._wp
            ! fd_coeff_s(1,  i) = -fd_coeff_s(-1,i)
            
            dx0 = s_cc(3) - s_cc(2)

            ! For the non‐uniform case i=1:  | 3/4*dx |  dx  |  dx  |
            ! (evaluation at the geometric center of the middle cell)
            fd_coeff_s(-1,1) = -64._wp/(105._wp*dx0)
            fd_coeff_s(0, 1) =  1._wp/(7._wp*dx0)
            fd_coeff_s(1, 1) =  7._wp/(15._wp*dx0)
            print *, 'DEBUG at 1', fd_coeff_s(-1:1,1), s_cc(2), s_cc(1), s_cc(0)
            
            ! For the non‐uniform case i=0:  | 1/2*dx | 3/4*dx |  dx  |
            ! (evaluation at the geometric center of the current cell, which here equals s_cc(0))
            fd_coeff_s(-1,0) = -14._wp/(15._wp*dx0)
            fd_coeff_s(0, 0) =  16._wp/(35._wp*dx0)
            fd_coeff_s(1, 0) =  10._wp/(21._wp*dx0)
            print *, 'DEBUG at 0', fd_coeff_s(-1:1,0), s_cc(1), s_cc(0), s_cc(-1)

            print *, lB, lE

            call exit(66)

            ! DEBUG at 1  -29866.666666666672        7000.0000000000009        22866.666666666672        4.0816326530612245E-005   2.0408163265306119E-005   5.1020408163265306E-006
            ! DEBUG at 0  -45733.333333333343        22400.000000000004        23333.333333333336        2.0408163265306119E-005   5.1020408163265306E-006  -5.1020408163265306E-006
            ! DEBUG at 1  -29866.666666666672        7000.0000000000009        22866.666666666672        4.0816326530612245E-005   2.0408163265306119E-005   5.1020408163265306E-006
            ! DEBUG at 0  -45733.333333333343        22400.000000000004        23333.333333333336        2.0408163265306119E-005   5.1020408163265306E-006  -5.1020408163265306E-006
        end if
        

    end subroutine s_compute_finite_difference_coefficients

end module m_finite_differences

