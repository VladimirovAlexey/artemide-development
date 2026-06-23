!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Test suite for aTMDe_Ogata
!
!   Tests the OgataIntegrator (Transform procedure) against known exact integrals.
!   Transform(n, F, q) computes  (1/2) int_0^infty b^(n+1) J_n(b*q) F(b) db
!   so all integrals below are retrieved as  2 * Transform(n, F_wrapped, q).
!
!   (1)  int b   J0(qb) exp(-ab)      = a / (a^2+q^2)^(3/2)
!        -> n=0, F = exp(-ab)
!
!   (2)  int b^2 J1(qb) exp(-ab)      = 3aq / (a^2+q^2)^(5/2)
!        -> n=1, F = exp(-ab)
!
!   (3)  int b^3 J2(qb) exp(-ab)      = 15aq^2 / (a^2+q^2)^(7/2)
!        -> n=2, F = exp(-ab)
!
!   (4)  int b   J0(qb) exp(-ab^2)    = exp(-q^2/4a) / (2a)
!        -> n=0, F = exp(-ab^2)
!
!   (5)  int b^2 J1(qb) exp(-ab^2)    = exp(-q^2/4a) * q / (4a^2)
!        -> n=1, F = exp(-ab^2)
!
!   (6)  int b^3 J2(qb) exp(-ab^2)    = exp(-q^2/4a) * q^2 / (8a^3)
!        -> n=2, F = exp(-ab^2)
!
!   (7)  int b^2 J0(qb) exp(-ab)      = (2a^2 - q^2) / (a^2+q^2)^(5/2)
!        -> n=0, F = b*exp(-ab)     [b^1 * F = b^2 * exp(-ab)]
!
!   Checked for a in {0.1, 1.0, 10.0}  and  q in {0.5,1,2,5,10,20,40}.
!   Output goes to terminal and Tests/aTMDe_Ogata/test.out.
!
!   Notes on Ogata quadrature behavior:
!   - hOGATA must be of order 10^-3 for these slowly-decaying integrands (exp(-ab));
!     larger h causes nodes to crowd near zero and misses the tail.
!   - At large q (q >> 1/a) the Bessel integrand oscillates faster than the
!     exponential damps it; convergence degrades and NC markers are expected.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_Ogata
use aTMDe_Ogata
use aTMDe_numerics
implicit none

!! local abstract interfaces so LoopAQ can accept any compatible procedure pair
abstract interface
    function integrand_t(b)
        import dp
        real(dp), intent(in) :: b
        real(dp) :: integrand_t
    end function integrand_t
    function exact_t(a, q)
        import dp
        real(dp), intent(in) :: a, q
        real(dp) :: exact_t
    end function exact_t
end interface

integer  :: fu
real(dp) :: a_par   !! shared integrand parameter — accessed by contained functions

real(dp), parameter :: avals(3) = [0.1_dp, 1._dp, 10._dp]
real(dp), parameter :: qvals(7) = [0.5_dp, 1._dp, 2._dp, 5._dp, 10._dp, 20._dp, 40._dp]

type(OgataIntegrator) :: ogata

open(newunit=fu, file='Tests/aTMDe_Ogata/test.out', status='replace', action='write')

!! single integrator instance; ww/bb tables are built for all Bessel orders 0..3
!! hOGATA ~ 10^-3 is required for exp(-ab) type integrands; larger h loses the tail
ogata = OgataIntegrator('test_Ogata', outputLevel=0, order_in=0, &
    tolerance_in=1.e-6_dp, hOGATA_in=0.001_dp, TMDmass_in=1._dp, kTmin=1.e-3_dp)

call SectionHeader('Lorentz decay  F(b) = exp(-a*b)')

call IntHeader('(1)  int b   J0(qb) exp(-ab)  =  a / (a^2+q^2)^(3/2)')
call LoopAQ(0, f_lor, ex1)

call IntHeader('(2)  int b^2 J1(qb) exp(-ab)  =  3aq / (a^2+q^2)^(5/2)')
call LoopAQ(1, f_lor, ex2)

call IntHeader('(3)  int b^3 J2(qb) exp(-ab)  =  15aq^2 / (a^2+q^2)^(7/2)')
call LoopAQ(2, f_lor, ex3)

call IntHeader('(7)  int b^2 J0(qb) exp(-ab)  =  (2a^2-q^2)/(a^2+q^2)^(5/2)   [F = b*exp(-ab)]')
call LoopAQ(0, f_lor_b, ex7)

call SectionHeader('Gaussian decay  F(b) = exp(-a*b^2)')

call IntHeader('(4)  int b   J0(qb) exp(-ab^2)  =  exp(-q^2/4a) / (2a)')
call LoopAQ(0, f_gauss, ex4)

call IntHeader('(5)  int b^2 J1(qb) exp(-ab^2)  =  exp(-q^2/4a) * q / (4a^2)')
call LoopAQ(1, f_gauss, ex5)

call IntHeader('(6)  int b^3 J2(qb) exp(-ab^2)  =  exp(-q^2/4a) * q^2 / (8a^3)   [F = exp(-ab^2)/b]')
call LoopAQ(2, f_gauss, ex6)

close(fu)

contains

!=================================================================================
!! integrand wrappers — all access a_par via host association

function f_lor(b)
    real(dp), intent(in) :: b
    real(dp) :: f_lor
    f_lor = exp(-a_par * b)
end function f_lor

function f_lor_b(b)               !! b*exp(-ab): gives int b^2 J0 via Transform n=0
    real(dp), intent(in) :: b
    real(dp) :: f_lor_b
    f_lor_b = b * exp(-a_par * b)
end function f_lor_b

function f_gauss(b)
    real(dp), intent(in) :: b
    real(dp) :: f_gauss
    f_gauss = exp(-a_par * b**2)
end function f_gauss

!=================================================================================
!! exact analytic formulas

function ex1(a, q)
    real(dp), intent(in) :: a, q; real(dp) :: ex1
    ex1 = a / (a**2 + q**2)**1.5_dp
end function ex1

function ex2(a, q)
    real(dp), intent(in) :: a, q; real(dp) :: ex2
    ex2 = 3._dp * a * q / (a**2 + q**2)**2.5_dp
end function ex2

function ex3(a, q)
    real(dp), intent(in) :: a, q; real(dp) :: ex3
    ex3 = 15._dp * a * q**2 / (a**2 + q**2)**3.5_dp
end function ex3

function ex4(a, q)
    real(dp), intent(in) :: a, q; real(dp) :: ex4
    ex4 = exp(-q**2 / (4._dp * a)) / (2._dp * a)
end function ex4

function ex5(a, q)
    real(dp), intent(in) :: a, q; real(dp) :: ex5
    ex5 = exp(-q**2 / (4._dp * a)) * q / (4._dp * a**2)
end function ex5

function ex6(a, q)
    real(dp), intent(in) :: a, q; real(dp) :: ex6
    ex6 = exp(-q**2 / (4._dp * a)) * q**2 / (8._dp * a**3)
end function ex6

function ex7(a, q)
    real(dp), intent(in) :: a, q; real(dp) :: ex7
    ex7 = (2._dp * a**2 - q**2) / (a**2 + q**2)**2.5_dp
end function ex7

!=================================================================================
!! test driver: loops over all (a, q) pairs

subroutine LoopAQ(n, f_func, ex_func)
    integer, intent(in) :: n
    procedure(integrand_t) :: f_func
    procedure(exact_t)     :: ex_func

    integer  :: ia, iq
    real(dp) :: q, computed, exact_val, ae, re
    logical  :: converged
    character(len=3) :: mark

    call ColHeader()
    do ia = 1, size(avals)
        a_par = avals(ia)
        do iq = 1, size(qvals)
            q = qvals(iq)
            call ogata%Transform(f_func, q, n, computed, converged)
            computed  = 2._dp * computed
            exact_val = ex_func(a_par, q)
            ae = abs(computed - exact_val)
            re = ae / max(abs(exact_val), 1.e-32_dp)
            mark = merge('   ', ' NC', converged)
            call PrintRow(a_par, q, computed, exact_val, ae, re, mark)
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
    call OutLine(repeat('=', 96))
    call OutLine('  ' // title)
    call OutLine(repeat('=', 96))
end subroutine SectionHeader

subroutine IntHeader(label)
    character(len=*), intent(in) :: label
    call OutLine('')
    call OutLine('  -- ' // label // ' --')
    call OutLine('')
end subroutine IntHeader

subroutine ColHeader()
    call OutLine('        a        q              Computed                  Exact' // &
                 '          |Error|    Rel.Err')
    call OutLine('   ' // repeat('-', 93))
end subroutine ColHeader

subroutine PrintRow(a, q, computed, exact, ae, re, mark)
    real(dp),         intent(in) :: a, q, computed, exact, ae, re
    character(len=3), intent(in) :: mark
    character(len=140) :: line
    write(line, '(3X,F7.2,2X,F7.2,2X,ES22.14,2X,ES22.14,2X,ES10.3,2X,ES10.3,A3)') &
        a, q, computed, exact, ae, re, mark
    call OutLine(trim(line))
end subroutine PrintRow

end program test_Ogata
