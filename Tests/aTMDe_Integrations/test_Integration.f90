!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Test program for aTMDe_Integration
!
!       Computes standard test integrals with known analytic values and prints
!       a comparison table: computed result, exact value, absolute and relative error.
!       Output is written to both stdout and test_Integration.out.
!
!       Test integrals (domain [0,1] unless stated):
!         f1: x^4              exact = 1/5
!         f2: exp(x)           exact = e - 1
!         f3: sin(pi*x)        exact = 2/pi
!         f4: 1/(1+x^2)        exact = pi/4
!         f5: cos(20*pi*x)     exact = 0  (oscillatory)
!
!       2D test integrals (domain [0,1]^2):
!         g1: x^2 + y^2            exact = 2/3
!         g2: sin(pi*x)*sin(pi*y)  exact = 4/pi^2
!         g3: exp(x+y)             exact = (e-1)^2
!
!       Array test (GK_array5):
!         h(x,k) = k*sin(pi*x),  exact integral = k * 2/pi,  k = -5..5
!
!                               A.Vladimirov (22.06.2026)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_Integration
use aTMDe_Numerics
use aTMDe_Integration
implicit none

character(len=*), parameter :: outfile = 'Tests/aTMDe_Integrations/test.out'
integer :: fu                          !! file unit, accessible in all contained procedures
real(dp), parameter :: tol = 1.0e-8_dp

real(dp) :: exact_poly, exact_exp, exact_sin, exact_rat, exact_osc
real(dp) :: exact_2d_quad, exact_2d_sinsin, exact_2d_exp
real(dp), dimension(-5:5) :: resArr
integer :: k

!!! analytic values
exact_poly      = 1.0_dp / 5.0_dp
exact_exp       = exp(1.0_dp) - 1.0_dp
exact_sin       = 2.0_dp / pi
exact_rat       = pi / 4.0_dp
exact_osc       = 0.0_dp
exact_2d_quad   = 2.0_dp / 3.0_dp
exact_2d_sinsin = 4.0_dp / pi2
exact_2d_exp    = (exp(1.0_dp) - 1.0_dp)**2

open(newunit=fu, file=outfile, status='replace', action='write')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  1D fixed-point methods
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call SectionHeader('1D fixed-point methods  [0,1]')

call MethodHeader('S5  (5-point composite Simpson)')
call Row('x^4',         Integrate_S5(f_poly, 0._dp, 1._dp), exact_poly)
call Row('exp(x)',       Integrate_S5(f_exp,  0._dp, 1._dp), exact_exp)
call Row('sin(pi*x)',    Integrate_S5(f_sin,  0._dp, 1._dp), exact_sin)
call Row('1/(1+x^2)',    Integrate_S5(f_rat,  0._dp, 1._dp), exact_rat)
call Row('cos(20pi*x)',  Integrate_S5(f_osc,  0._dp, 1._dp), exact_osc)

call MethodHeader('G3  (3-point Gauss,  exact for poly to deg 5)')
call Row('x^4',         Integrate_G3(f_poly, 0._dp, 1._dp), exact_poly)
call Row('exp(x)',       Integrate_G3(f_exp,  0._dp, 1._dp), exact_exp)
call Row('sin(pi*x)',    Integrate_G3(f_sin,  0._dp, 1._dp), exact_sin)
call Row('1/(1+x^2)',    Integrate_G3(f_rat,  0._dp, 1._dp), exact_rat)
call Row('cos(20pi*x)',  Integrate_G3(f_osc,  0._dp, 1._dp), exact_osc)

call MethodHeader('G7  (7-point Gauss,  exact for poly to deg 13)')
call Row('x^4',         Integrate_G7(f_poly, 0._dp, 1._dp), exact_poly)
call Row('exp(x)',       Integrate_G7(f_exp,  0._dp, 1._dp), exact_exp)
call Row('sin(pi*x)',    Integrate_G7(f_sin,  0._dp, 1._dp), exact_sin)
call Row('1/(1+x^2)',    Integrate_G7(f_rat,  0._dp, 1._dp), exact_rat)
call Row('cos(20pi*x)',  Integrate_G7(f_osc,  0._dp, 1._dp), exact_osc)

call MethodHeader('K15  (15-point Kronrod,  exact for poly to deg 29)')
call Row('x^4',         Integrate_K15(f_poly, 0._dp, 1._dp), exact_poly)
call Row('exp(x)',       Integrate_K15(f_exp,  0._dp, 1._dp), exact_exp)
call Row('sin(pi*x)',    Integrate_K15(f_sin,  0._dp, 1._dp), exact_sin)
call Row('1/(1+x^2)',    Integrate_K15(f_rat,  0._dp, 1._dp), exact_rat)
call Row('cos(20pi*x)',  Integrate_K15(f_osc,  0._dp, 1._dp), exact_osc)

call MethodHeader('K21  (21-point Kronrod,  exact for poly to deg 41)')
call Row('x^4',         Integrate_K21(f_poly, 0._dp, 1._dp), exact_poly)
call Row('exp(x)',       Integrate_K21(f_exp,  0._dp, 1._dp), exact_exp)
call Row('sin(pi*x)',    Integrate_K21(f_sin,  0._dp, 1._dp), exact_sin)
call Row('1/(1+x^2)',    Integrate_K21(f_rat,  0._dp, 1._dp), exact_rat)
call Row('cos(20pi*x)',  Integrate_K21(f_osc,  0._dp, 1._dp), exact_osc)

call MethodHeader('K41  (41-point Kronrod,  exact for poly to deg 81)')
call Row('x^4',         Integrate_K41(f_poly, 0._dp, 1._dp), exact_poly)
call Row('exp(x)',       Integrate_K41(f_exp,  0._dp, 1._dp), exact_exp)
call Row('sin(pi*x)',    Integrate_K41(f_sin,  0._dp, 1._dp), exact_sin)
call Row('1/(1+x^2)',    Integrate_K41(f_rat,  0._dp, 1._dp), exact_rat)
call Row('cos(20pi*x)',  Integrate_K41(f_osc,  0._dp, 1._dp), exact_osc)

call MethodHeader('SN  (N-point composite Simpson,  convergence on exp(x))')
call Row('N=10',   Integrate_SN(f_exp, 0._dp, 1._dp, 10),   exact_exp)
call Row('N=100',  Integrate_SN(f_exp, 0._dp, 1._dp, 100),  exact_exp)
call Row('N=1000', Integrate_SN(f_exp, 0._dp, 1._dp, 1000), exact_exp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  1D adaptive methods
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call SectionHeader('1D adaptive methods  [0,1],  tolerance = 1E-8')

call MethodHeader('SA  (adaptive Simpson)')
call Row('x^4',        Integrate_SA(f_poly, 0._dp, 1._dp, tol), exact_poly)
call Row('exp(x)',      Integrate_SA(f_exp,  0._dp, 1._dp, tol), exact_exp)
call Row('sin(pi*x)',   Integrate_SA(f_sin,  0._dp, 1._dp, tol), exact_sin)
call Row('1/(1+x^2)',   Integrate_SA(f_rat,  0._dp, 1._dp, tol), exact_rat)
call Row('cos(20pi*x)', Integrate_SA(f_osc,  0._dp, 1._dp, tol), exact_osc)

call MethodHeader('GK  (adaptive Gauss-Kronrod 7/15)')
call Row('x^4',        Integrate_GK(f_poly, 0._dp, 1._dp, tol), exact_poly)
call Row('exp(x)',      Integrate_GK(f_exp,  0._dp, 1._dp, tol), exact_exp)
call Row('sin(pi*x)',   Integrate_GK(f_sin,  0._dp, 1._dp, tol), exact_sin)
call Row('1/(1+x^2)',   Integrate_GK(f_rat,  0._dp, 1._dp, tol), exact_rat)
call Row('cos(20pi*x)', Integrate_GK(f_osc,  0._dp, 1._dp, tol), exact_osc)

call MethodHeader('GK2041  (adaptive Gauss-Kronrod 20/41)')
call Row('x^4',        Integrate_GK2041(f_poly, 0._dp, 1._dp, tol), exact_poly)
call Row('exp(x)',      Integrate_GK2041(f_exp,  0._dp, 1._dp, tol), exact_exp)
call Row('sin(pi*x)',   Integrate_GK2041(f_sin,  0._dp, 1._dp, tol), exact_sin)
call Row('1/(1+x^2)',   Integrate_GK2041(f_rat,  0._dp, 1._dp, tol), exact_rat)
call Row('cos(20pi*x)', Integrate_GK2041(f_osc,  0._dp, 1._dp, tol), exact_osc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  GK array (-5:5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call SectionHeader('GK_array5  [0,1],  f(x,k) = k*sin(pi*x),  exact = k*2/pi,  tol = 1E-8')
call Out('(3X,A6,2X,A20,2X,A20,2X,A12,2X,A10)', &
         'Index','Computed','Exact','|Error|','Rel.Err')
call OutLine(repeat('-',75))

resArr = Integrate_GK_array5(f_arr, 0._dp, 1._dp, tol)
do k = -5, 5
    call RowArr(k, resArr(k), real(k,dp)*exact_sin)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  2D methods
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call SectionHeader('2D methods  [0,1]^2,  tolerance = 1E-8')

call MethodHeader('SA_2D  (adaptive 2D Simpson)')
call Row('x^2+y^2',       Integrate_SA_2D(f2d_quad,   0._dp,1._dp, 0._dp,1._dp, tol), exact_2d_quad)
call Row('sin(px)sin(py)', Integrate_SA_2D(f2d_sinsin, 0._dp,1._dp, 0._dp,1._dp, tol), exact_2d_sinsin)
call Row('exp(x+y)',       Integrate_SA_2D(f2d_exp,    0._dp,1._dp, 0._dp,1._dp, tol), exact_2d_exp)

call MethodHeader('Stroud35_55  (adaptive 2D Stroud 3-5 / 5-5)')
call Row('x^2+y^2',       Integrate2D_Stroud35_55(f2d_quad,   0._dp,1._dp, 0._dp,1._dp, tol), exact_2d_quad)
call Row('sin(px)sin(py)', Integrate2D_Stroud35_55(f2d_sinsin, 0._dp,1._dp, 0._dp,1._dp, tol), exact_2d_sinsin)
call Row('exp(x+y)',       Integrate2D_Stroud35_55(f2d_exp,    0._dp,1._dp, 0._dp,1._dp, tol), exact_2d_exp)

call MethodHeader('Stroud35_56  (adaptive 2D Stroud 3-5 / 5-6)')
call Row('x^2+y^2',       Integrate2D_Stroud35_56(f2d_quad,   0._dp,1._dp, 0._dp,1._dp, tol), exact_2d_quad)
call Row('sin(px)sin(py)', Integrate2D_Stroud35_56(f2d_sinsin, 0._dp,1._dp, 0._dp,1._dp, tol), exact_2d_sinsin)
call Row('exp(x+y)',       Integrate2D_Stroud35_56(f2d_exp,    0._dp,1._dp, 0._dp,1._dp, tol), exact_2d_exp)

call OutLine(repeat('=',92))
call OutLine('')

close(fu)
write(*,'(A)') '  Output also written to '//outfile

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  1D test functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function f_poly(x) result(r)
    real(dp), intent(in) :: x
    real(dp) :: r
    r = x**4
end function

function f_exp(x) result(r)
    real(dp), intent(in) :: x
    real(dp) :: r
    r = exp(x)
end function

function f_sin(x) result(r)
    real(dp), intent(in) :: x
    real(dp) :: r
    r = sin(pi*x)
end function

function f_rat(x) result(r)
    real(dp), intent(in) :: x
    real(dp) :: r
    r = 1.0_dp / (1.0_dp + x**2)
end function

function f_osc(x) result(r)
    real(dp), intent(in) :: x
    real(dp) :: r
    r = cos(20.0_dp*pi*x)
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Array test function  h(x,k) = k*sin(pi*x)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function f_arr(x) result(r)
    real(dp), intent(in) :: x
    real(dp), dimension(-5:5) :: r
    integer :: j
    do j = -5, 5
        r(j) = real(j,dp) * sin(pi*x)
    end do
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  2D test functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function f2d_quad(x,y) result(r)
    real(dp), intent(in) :: x, y
    real(dp) :: r
    r = x**2 + y**2
end function

function f2d_sinsin(x,y) result(r)
    real(dp), intent(in) :: x, y
    real(dp) :: r
    r = sin(pi*x) * sin(pi*y)
end function

function f2d_exp(x,y) result(r)
    real(dp), intent(in) :: x, y
    real(dp) :: r
    r = exp(x + y)
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Output helpers  (write to stdout and file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine OutLine(str)
    character(len=*), intent(in) :: str
    write(*,'(A)') str
    write(fu,'(A)') str
end subroutine

subroutine Out(fmt, a1, a2, a3, a4, a5)
    character(len=*), intent(in) :: fmt, a1, a2, a3, a4, a5
    write(*,fmt) a1, a2, a3, a4, a5
    write(fu,fmt) a1, a2, a3, a4, a5
end subroutine

subroutine SectionHeader(title)
    character(len=*), intent(in) :: title
    call OutLine('')
    call OutLine(repeat('=',92))
    call OutLine('  '//title)
    call OutLine(repeat('=',92))
    write(*,'(3X,A20,2X,A20,2X,A20,2X,A12,2X,A10)') &
        'Function','Computed','Exact','|Error|','Rel.Err'
    write(fu,'(3X,A20,2X,A20,2X,A20,2X,A12,2X,A10)') &
        'Function','Computed','Exact','|Error|','Rel.Err'
    call OutLine('   '//repeat('-',88))
end subroutine

subroutine MethodHeader(title)
    character(len=*), intent(in) :: title
    call OutLine('')
    call OutLine('  -- '//trim(title)//' --')
end subroutine

subroutine Row(fname, computed, exact)
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: computed, exact
    real(dp) :: ae, re
    ae = abs(computed - exact)
    if(abs(exact) > epsilon(1.0_dp)) then
        re = ae / abs(exact)
    else
        re = ae
    end if
    write(*,'(3X,A20,2X,ES20.12,2X,ES20.12,2X,ES12.4,2X,ES10.3)') &
        fname, computed, exact, ae, re
    write(fu,'(3X,A20,2X,ES20.12,2X,ES20.12,2X,ES12.4,2X,ES10.3)') &
        fname, computed, exact, ae, re
end subroutine

subroutine RowArr(idx, computed, exact)
    integer, intent(in) :: idx
    real(dp), intent(in) :: computed, exact
    real(dp) :: ae, re
    character(len=6) :: s
    ae = abs(computed - exact)
    if(abs(exact) > epsilon(1.0_dp)) then
        re = ae / abs(exact)
    else
        re = ae
    end if
    write(s,'(I4)') idx
    write(*,'(3X,A6,2X,ES20.12,2X,ES20.12,2X,ES12.4,2X,ES10.3)') &
        adjustl(s), computed, exact, ae, re
    write(fu,'(3X,A6,2X,ES20.12,2X,ES20.12,2X,ES12.4,2X,ES10.3)') &
        adjustl(s), computed, exact, ae, re
end subroutine

end program test_Integration
