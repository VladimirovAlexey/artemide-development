!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Test suite for aTMDe_optGrid
!
!   Builds an optGrid using uTMDPDF_OPE_convolution(x,b,h,.false.) as the test function,
!   then compares extraction from the grid against direct OPE computation.
!
!   Error metric: ae / f_scale(x,fv)  where f_scale = max_b |f(x,b,fv)|
!   This normalises by the function scale at each x, so errors at deep minima
!   (where |f| << f_scale) are not artificially inflated in the relative sense.
!
!   Reported accuracy: "accurate digits" = -log10( ae / f_scale )
!   Negative value means even the order of magnitude is wrong.
!
!   The ini-file must have grid preparation disabled (*p2 = F in the *B section)
!   so that uTMDPDF_OPE_convolution calls CxF_compute directly (the "exact" reference).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_optGrid
use aTMDe_IO
use aTMDe_optGrid
use aTMDe_numerics
use aTMDe_control
use uTMDPDF_OPE
implicit none

!! ---- configuration (update for your setup) ----
character(*), parameter :: inifile   = 'test.atmde'
character(*), parameter :: prefix    = 'Tests/aTMDe_optGrid/'
real(dp),     parameter :: threshold = 1.e-5_dp   !! flag when ae/f_scale exceeds this

!! ---- test grid: points chosen to fall between Chebyshev nodes ----
integer,  parameter :: Nx = 45, Nb = 45
real(dp), parameter :: xMIN = 0.00001_dp
real(dp), parameter :: bMAX = 15._dp
real(dp) :: xtest(Nx), btest(Nb)

!! ---- number of timing repetitions ----
integer, parameter :: Nrep = 20

!! ---- pre-computed exact values and their per-x scale ----
real(dp) :: exact_store(Nx, Nb, -5:5)
real(dp) :: f_scale(Nx, -5:5)       !! max_b |f(x,b,fv)| for each (x,fv)

integer  :: fu, ix, ib, fv, n_bad, irep
real(dp) :: time1, time2, time_make, time_extract, time_direct
real(dp), dimension(-5:5) :: from_grid, from_exact
real(dp) :: ae, re, nd, nd_min
character(len=200) :: buf

type(optGrid) :: grid

do ix = 1, Nx
    xtest(ix) = 10._dp**(log10(xMIN) * (Nx+1-ix) / (Nx+1))
end do
do ib = 1, Nb
    btest(ib) = bMAX * ib / (Nb+1)
end do

open(newunit=fu, file='Tests/aTMDe_optGrid/test.out', status='replace', action='write')

!! ============================================================
!! Initialization
!! ============================================================
call artemide_Initialize(inifile, prefix=prefix)
call artemide_SetNPparameters_uTMDPDF([&
    0.874245_dp, 0.913883_dp, 0.991563_dp, 6.05412_dp, &
    0.353908_dp, 46.6064_dp,  0.115161_dp, 1.53235_dp, &
    1.31966_dp,  0.434833_dp, 0.0_dp,      0.0_dp])

!! ============================================================
!! Build grid and measure construction time
!! ============================================================
call cpu_time(time1)
grid = optGrid(trim(prefix)//trim(inifile), '*4   ', '*E   ', 1, .false., 'test_optGrid', 1)
call grid%MakeGrid(F_test)
call cpu_time(time2)
time_make = time2 - time1

!! ============================================================
!! Pre-compute all exact values; derive per-x scale
!! ============================================================
f_scale = 0._dp
do ix = 1, Nx
    do ib = 1, Nb
        from_exact = x_times_OPE(xtest(ix), btest(ib), 1)
        exact_store(ix, ib, -5:5) = from_exact
        do fv = -5, 5
            if (fv == 0) cycle
            f_scale(ix, fv) = max(f_scale(ix, fv), abs(from_exact(fv)))
        end do
    end do
end do

!! ============================================================
!! Accuracy test
!! ============================================================
call PutLine(repeat('=', 95))
call PutLine('  optGrid accuracy test:  grid extraction vs. direct OPE computation')
call PutLine('  Error metric: ae / max_b|f(x,fv)|   (normalised by function scale at each x)')
call PutLine('  Flagged when normalised error > '//trim(numToStr(threshold)))
call PutLine(repeat('=', 95))
call PutLine('')
call PutLine('  x            b          fl   from_grid             from_exact' &
             //'            |error|    acc.dig')
call PutLine('  '//repeat('-', 98))

n_bad  = 0
nd_min = 1.e30_dp

do ix = 1, Nx
    do ib = 1, Nb
        from_grid = grid%Extract(xtest(ix), btest(ib), 1)
        do fv = -5, 5
            if (fv == 0) cycle
            ae = abs(from_grid(fv) - exact_store(ix, ib, fv))
            re = ae / max(f_scale(ix, fv), 1.e-32_dp)
            nd = -log10(max(re, 1.e-16_dp))
            if (nd < nd_min) nd_min = nd
            if (re > threshold) then
                n_bad = n_bad + 1
                write(buf, '(2X,ES10.3,2X,ES10.3,2X,I3,2X,ES22.14,2X,ES22.14,2X,ES10.3,2X,F7.2)') &
                    xtest(ix), btest(ib), fv, from_grid(fv), exact_store(ix,ib,fv), ae, nd
                call PutLine(trim(buf))
            end if
        end do
    end do
end do

call PutLine('')

!! ============================================================
!! Timing: pure extraction (Nrep * Nx * Nb calls)
!! ============================================================
call cpu_time(time1)
do irep = 1, Nrep
    do ix = 1, Nx
        do ib = 1, Nb
            from_grid = grid%Extract(xtest(ix), btest(ib), 1)
        end do
    end do
end do
call cpu_time(time2)
time_extract = time2 - time1

!! Timing: direct OPE computation (same number of calls)
call cpu_time(time1)
do irep = 1, Nrep
    do ix = 1, Nx
        do ib = 1, Nb
            from_exact = x_times_OPE(xtest(ix), btest(ib), 1)
        end do
    end do
end do
call cpu_time(time2)
time_direct = time2 - time1

!! ============================================================
!! Summary
!! ============================================================
call PutLine(repeat('=', 95))
call PutLine('  Summary')
call PutLine(repeat('=', 95))
write(buf, '("  Test points (x x b):             ",I3," x ",I3," = ",I5)') Nx, Nb, Nx*Nb
call PutLine(trim(buf))
write(buf, '("  Flavors tested:                  -5..5  (10 per point, total ",I6,")")') Nx*Nb*10
call PutLine(trim(buf))
write(buf, '("  Flagged cases (norm. err > ",ES8.1,"):  ",I5)') threshold, n_bad
call PutLine(trim(buf))
write(buf, '("  Worst accurate digits:           ",F7.2," (= min over all points)")') nd_min
call PutLine(trim(buf))
call PutLine('')
write(buf, '("  MakeGrid time:        ",F10.4," s")') time_make
call PutLine(trim(buf))
write(buf, '("  Grid extraction:      ",F10.4," s  (",I6," calls)")') &
    time_extract, Nrep*Nx*Nb
call PutLine(trim(buf))
write(buf, '("  Direct OPE:           ",F10.4," s  (",I6," calls)")') &
    time_direct, Nrep*Nx*Nb
call PutLine(trim(buf))
if (time_extract > 0._dp) then
    write(buf, '("  Speedup factor:       ",F10.1,"x")') time_direct / time_extract
    call PutLine(trim(buf))
end if

close(fu)

contains

function F_test(x, b, h)
    real(dp), intent(in) :: x, b
    integer,  intent(in) :: h
    real(dp), dimension(-5:5) :: F_test
    F_test = x * uTMDPDF_OPE_convolution(x, b, h, .false.)
end function F_test

function x_times_OPE(x, b, h)
    real(dp), intent(in) :: x, b
    integer,  intent(in) :: h
    real(dp), dimension(-5:5) :: x_times_OPE
    x_times_OPE = x * uTMDPDF_OPE_convolution(x, b, h, .false.)
end function x_times_OPE

subroutine PutLine(str)
    character(len=*), intent(in) :: str
    write(*,  '(A)') str
    write(fu, '(A)') str
end subroutine PutLine

end program test_optGrid
