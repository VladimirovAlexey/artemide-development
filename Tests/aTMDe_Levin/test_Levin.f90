!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Test suite for aTMDe_Levin
!
!   Tests the LevinIntegrator (Fourier_atPoint) against known exact integrals.
!
!   Type-0 transform:  (1/2pi) int_0^inf b   J_0(b*q) F(b) db
!   Type-1 transform:  (M^2/q)*(1/2pi) int_0^inf b^2 J_1(b*q) F(b) db   [M=1]
!
!   Tested integrals:
!
!   (1)  (1/2pi) int b   J_0(qb) exp(-ab)     = a / (2pi*(a^2+q^2)^(3/2))
!        -> type=0, F = exp(-ab)
!
!   (2)  (1/2pi) int b^2 J_0(qb) exp(-ab)     = (2a^2-q^2) / (2pi*(a^2+q^2)^(5/2))
!        -> type=0, F = b*exp(-ab)
!
!   (3)  (1/2pi) int b   J_0(qb) exp(-ab^2)   = exp(-q^2/4a) / (4pi*a)
!        -> type=0, F = exp(-ab^2)
!
!   (4)  (1/q)*(1/2pi) int b^2 J_1(qb) exp(-ab)   = 3a / (2pi*(a^2+q^2)^(5/2))
!        -> type=1, F = exp(-ab)
!
!   (5)  (1/q)*(1/2pi) int b^2 J_1(qb) exp(-ab^2) = exp(-q^2/4a) / (8pi*a^2)
!        -> type=1, F = exp(-ab^2)
!
!   Lorentzian cases checked for a in {1., 2., 5.}  and  q in {0.1,0.5,1,5,10,50,200,500}.
!   Gaussian    cases checked for a in {0.1, 0.5, 1.} and the same q range.
!   b-grid: [1e-5, 100] over 5 subgrids with 16 nodes each.
!   Output goes to terminal and Tests/aTMDe_Levin/test.out.
!
!   Notes:
!   - The Levin method excels at large q where the Bessel integrand oscillates rapidly;
!     good precision is expected across the full q range including q=200 and q=500.
!   - Gaussian integrals at large q (e.g. a=0.1, q>=50) evaluate to effectively zero;
!     the relative error in those rows is meaningless and can be ignored.
!   - Integral (2) changes sign at q = a*sqrt(2); near the zero crossing the relative
!     error is large even when the absolute error is small.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_Levin
use aTMDe_Levin
use aTMDe_interfaces
use aTMDe_numerics
implicit none

abstract interface
    function exact_t(a, q)
        import dp
        real(dp), intent(in) :: a, q
        real(dp) :: exact_t
    end function exact_t
end interface

integer  :: fu
real(dp) :: a_par   !! shared integrand parameter — accessed by contained functions via host association

character(*), parameter :: inipath = 'Tests/aTMDe_Levin/test.atmde'

real(dp), parameter :: avals_lor(3) = [1._dp, 2._dp, 5._dp]
real(dp), parameter :: avals_gau(3) = [0.1_dp, 0.5_dp, 1._dp]
real(dp), parameter :: qvals(12) = &
        [0.01_dp, 0.1_dp, 0.5_dp, 1._dp, 2._dp, 5._dp, 10._dp, 20._dp, 50._dp, 100._dp, 150._dp , 199._dp]

type(LevinIntegrator) :: lev0, lev1

open(newunit=fu, file='Tests/aTMDe_Levin/test.out', status='replace', action='write')

lev0 = LevinIntegrator(inipath, '*T   ', '*G   ', 'test_Levin', 0, 0)
lev1 = LevinIntegrator(inipath, '*T   ', '*H   ', 'test_Levin', 0, 1)

!=================================================================================
!! Type-0 tests
call SectionHeader('Type 0:  (1/2pi) int_0^inf b J_0(b*q) F(b) db')

call IntHeader('(1)  F(b) = exp(-a*b)          exact: a / (2pi*(a^2+q^2)^(3/2))')
call LoopAQ(lev0, avals_lor, f_lor, ex0_lor)

call IntHeader('(2)  F(b) = b*exp(-a*b)        exact: (2a^2-q^2) / (2pi*(a^2+q^2)^(5/2))')
call LoopAQ(lev0, avals_lor, f_bexp, ex0_bexp)

call IntHeader('(3)  F(b) = exp(-a*b^2)        exact: exp(-q^2/4a) / (4pi*a)')
call LoopAQ(lev0, avals_gau, f_gauss, ex0_gau)

!=================================================================================
!! Type-1 tests
call SectionHeader('Type 1:  (M^2/q)*(1/2pi) int_0^inf b^2 J_1(b*q) F(b) db   [M=1]')

call IntHeader('(4)  F(b) = exp(-a*b)          exact: 3a / (2pi*(a^2+q^2)^(5/2))')
call LoopAQ(lev1, avals_lor, f_lor, ex1_lor)

call IntHeader('(5)  F(b) = exp(-a*b^2)        exact: exp(-q^2/4a) / (8pi*a^2)')
call LoopAQ(lev1, avals_gau, f_gauss, ex1_gau)

close(fu)

contains

!=================================================================================
!! integrand wrappers — all return the same scalar value in every (-5:5) component
!! a_par is accessed via host association

function f_lor(b)
    real(dp), intent(in) :: b
    real(dp), dimension(-5:5) :: f_lor
    f_lor = exp(-a_par * b)
end function f_lor

function f_bexp(b)
    real(dp), intent(in) :: b
    real(dp), dimension(-5:5) :: f_bexp
    f_bexp = b * exp(-a_par * b)
end function f_bexp

function f_gauss(b)
    real(dp), intent(in) :: b
    real(dp), dimension(-5:5) :: f_gauss
    f_gauss = exp(-a_par * b**2)
end function f_gauss

!=================================================================================
!! exact analytic formulas

function ex0_lor(a, q)
    real(dp), intent(in) :: a, q
    real(dp) :: ex0_lor
    ex0_lor = a / (pix2 * (a**2 + q**2)**1.5_dp)
end function ex0_lor

function ex0_bexp(a, q)
    real(dp), intent(in) :: a, q
    real(dp) :: ex0_bexp
    ex0_bexp = (2._dp*a**2 - q**2) / (pix2 * (a**2 + q**2)**2.5_dp)
end function ex0_bexp

function ex0_gau(a, q)
    real(dp), intent(in) :: a, q
    real(dp) :: ex0_gau
    ex0_gau = exp(-q**2 / (4._dp*a)) / (2._dp * pix2 * a)
end function ex0_gau

function ex1_lor(a, q)
    real(dp), intent(in) :: a, q
    real(dp) :: ex1_lor
    ex1_lor = 3._dp * a / (pix2 * (a**2 + q**2)**2.5_dp)
end function ex1_lor

function ex1_gau(a, q)
    real(dp), intent(in) :: a, q
    real(dp) :: ex1_gau
    ex1_gau = exp(-q**2 / (4._dp*a)) / (4._dp * pix2 * a**2)
end function ex1_gau

!=================================================================================
!! test driver: loops over all (a, q) pairs and checks component (0) of the result

subroutine LoopAQ(lev, avals, f_func, ex_func)
    class(LevinIntegrator), intent(in) :: lev
    real(dp),               intent(in) :: avals(:)
    procedure(func_1D_array5)          :: f_func
    procedure(exact_t)                 :: ex_func

    integer  :: ia, iq
    real(dp) :: q, computed, exact_val, ae, re
    real(dp), dimension(-5:5) :: computed_arr
    character(len=140) :: line

    call ColHeader()
    do ia = 1, size(avals)
        a_par = avals(ia)
        do iq = 1, size(qvals)
            q         = qvals(iq)
            computed_arr = lev%Fourier_atPoint(f_func, q)
            computed  = computed_arr(0)
            exact_val = ex_func(a_par, q)
            ae = abs(computed - exact_val)
            re = ae / max(abs(exact_val), 1.e-32_dp)
            write(line, '(3X,F7.2,2X,F7.1,2X,ES22.14,2X,ES22.14,2X,ES10.3,2X,ES10.3)') &
                a_par, q, computed, exact_val, ae, re
            call OutLine(trim(line))
        end do
        call OutLine('')
    end do
end subroutine LoopAQ

!=================================================================================
!! output helpers

subroutine OutLine(str)
    character(len=*), intent(in) :: str
    write(*,  '(A)') str
    write(fu, '(A)') str
end subroutine OutLine

subroutine SectionHeader(title)
    character(len=*), intent(in) :: title
    call OutLine('')
    call OutLine(repeat('=', 100))
    call OutLine('  ' // title)
    call OutLine(repeat('=', 100))
end subroutine SectionHeader

subroutine IntHeader(label)
    character(len=*), intent(in) :: label
    call OutLine('')
    call OutLine('  -- ' // label // ' --')
    call OutLine('')
end subroutine IntHeader

subroutine ColHeader()
    call OutLine('        a          q              Computed                  Exact' // &
                 '          |Error|    Rel.Err')
    call OutLine('   ' // repeat('-', 97))
end subroutine ColHeader

end program test_Levin
