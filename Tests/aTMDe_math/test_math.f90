!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Test suite for aTMDe_math
!
!   Test 1: ChebyshevT_array via the exact identity T_n(cos theta) = cos(n theta)
!           Also verifies symmetry T_{-n} = T_n
!   Test 2: ChebyshevT_int_array vs. numerical definite integral via Integrate_GK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_math
use aTMDe_math
use aTMDe_Integration
use aTMDe_numerics
implicit none

integer :: fu
integer :: k_int   !! shared index for the integrand function

open(newunit=fu, file='Tests/aTMDe_math/test.out', status='replace', action='write')

call TestChebyshevT()
call TestChebyshevT_int()

close(fu)

contains

!=================================================================================
subroutine TestChebyshevT()
    integer,  parameter :: nvals(5)  = [5, 10, 20, 50, 100]
    real(dp), parameter :: thetas(4) = [pi/7._dp, pi/4._dp, pi/3._dp, 2._dp*pi/5._dp]
    integer  :: in, it, n
    real(dp) :: theta, x, computed, exact, ae, re
    real(dp), allocatable :: Tv(:)

    call SectionHeader('ChebyshevT_array:  T_n(cos theta) = cos(n*theta)')
    call OutLine('      n   theta/pi              Computed                  Exact          |Error|    Rel.Err')
    call OutLine('   '//repeat('-',93))

    do in = 1, size(nvals)
        n = nvals(in)
        allocate(Tv(0:n))
        do it = 1, size(thetas)
            theta = thetas(it)
            x     = cos(theta)
            Tv    = ChebyshevT_array(n, x)
            computed = Tv(n)
            exact    = cos(real(n,dp)*theta)
            ae = abs(computed - exact)
            re = ae / max(abs(exact), epsilon(1._dp))
            call RowCheby(n, theta/pi, computed, exact, ae, re)
        end do
        deallocate(Tv)
        call OutLine('')
    end do

    !! T_{-n}(x) must equal T_n(x)
    call OutLine('  -- Symmetry: T_{-n} = T_n --')
    call OutLine('      n   theta/pi                T_n                 T_{-n}             |Diff|')
    call OutLine('   '//repeat('-',82))
    do it = 1, size(thetas)
        theta = thetas(it)
        x = cos(theta)
        call RowSymm(20, theta/pi, ChebyshevT(20,x), ChebyshevT(-20,x))
    end do
    call OutLine('')

end subroutine TestChebyshevT

!=================================================================================
subroutine TestChebyshevT_int()
    integer,  parameter :: n   = 15
    real(dp), parameter :: tol = 1.e-10_dp
    real(dp), parameter :: xvals(3) = [0.3_dp, 0.7_dp, -0.5_dp]
    integer  :: ix, k
    real(dp) :: x, int_formula, int_numerical, ae, re
    real(dp) :: Fint_x(0:n), Fint_0(0:n)

    call SectionHeader('ChebyshevT_int_array:  integral of T_k from 0 to x  (tol=1E-10 GK)')
    call OutLine('      k              Formula               GK numerical        |Error|    Rel.Err')
    call OutLine('   '//repeat('-',80))

    do ix = 1, size(xvals)
        x = xvals(ix)
        Fint_x = ChebyshevT_int_array(n, x)
        Fint_0 = ChebyshevT_int_array(n, 0._dp)

        write(*,  '(3X,A,F5.2)') 'x = ', x
        write(fu, '(3X,A,F5.2)') 'x = ', x

        do k = 0, n
            k_int         = k
            int_formula   = Fint_x(k) - Fint_0(k)
            int_numerical = Integrate_GK(f_Tk, 0._dp, x, tol)
            ae = abs(int_formula - int_numerical)
            re = ae / max(abs(int_numerical), epsilon(1._dp))
            call RowInt(k, int_formula, int_numerical, ae, re)
        end do
        call OutLine('')
    end do

end subroutine TestChebyshevT_int

!=================================================================================
!! Integrand: T_{k_int}(t), uses host-associated k_int
function f_Tk(t)
    real(dp), intent(in) :: t
    real(dp) :: f_Tk
    real(dp) :: Tv(0:k_int)
    Tv   = ChebyshevT_array(k_int, t)
    f_Tk = Tv(k_int)
end function f_Tk

!=================================================================================
!! Output helpers

subroutine OutLine(str)
    character(len=*), intent(in) :: str
    write(*,  '(A)') str
    write(fu, '(A)') str
end subroutine OutLine

subroutine SectionHeader(title)
    character(len=*), intent(in) :: title
    call OutLine('')
    call OutLine(repeat('=',96))
    call OutLine('  '//title)
    call OutLine(repeat('=',96))
end subroutine SectionHeader

subroutine RowCheby(n, theta_pi, computed, exact, ae, re)
    integer,  intent(in) :: n
    real(dp), intent(in) :: theta_pi, computed, exact, ae, re
    character(len=130) :: line
    write(line,'(3X,I6,2X,F8.5,2X,ES22.14,2X,ES22.14,2X,ES10.3,2X,ES10.3)') &
        n, theta_pi, computed, exact, ae, re
    call OutLine(trim(line))
end subroutine RowCheby

subroutine RowSymm(n, theta_pi, Tpos, Tneg)
    integer,  intent(in) :: n
    real(dp), intent(in) :: theta_pi, Tpos, Tneg
    character(len=130) :: line
    write(line,'(3X,I6,2X,F8.5,2X,ES22.14,2X,ES22.14,2X,ES10.3)') &
        n, theta_pi, Tpos, Tneg, abs(Tpos-Tneg)
    call OutLine(trim(line))
end subroutine RowSymm

subroutine RowInt(k, formula, numerical, ae, re)
    integer,  intent(in) :: k
    real(dp), intent(in) :: formula, numerical, ae, re
    character(len=130) :: line
    write(line,'(3X,I6,4X,ES22.14,2X,ES22.14,2X,ES10.3,2X,ES10.3)') &
        k, formula, numerical, ae, re
    call OutLine(trim(line))
end subroutine RowInt

end program test_math
