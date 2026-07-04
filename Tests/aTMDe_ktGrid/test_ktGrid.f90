!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Test suite for aTMDe_ktGrid
!
!   Builds a ktGrid using uTMDPDF as the reference TMD module:
!     - toGrid_test fills the grid via Fourier_Levin_array(toFourier_inner)
!     - toFourier_inner evaluates uTMDPDF_inB(x, b, Q, Q^2, h) — the full evolved b-space TMD
!
!   This is exactly the integrand that uTMDPDF_inKT = TMD_ev_inKT integrates via
!   Fourier_Levin(toFourier, kT), so the two computations share the same mathematical
!   object and the test isolates interpolation error from physical modelling choices.
!
!   Reference: uTMDPDF_inKT(x, kT, Q, Q^2, h) = TMD_ev_inKT (direct Fourier_Levin)
!
!   Test kT points: log scale from kMIN to Q  (kMIN = kRanges(0) from the config *F / *p7)
!   Test Q values:  2, 5, 10, 50, 100  GeV
!   Test x points:  log scale from xMIN_test to ~0.9
!
!   Error metric: ae / max_kT|f(x,fv)|  (normalised by function scale per (x,Q,fv))
!   Reported accuracy: "accurate digits" = -log10( ae / f_scale )
!
!   Config requirement: *p1 = T in the *F section so that
!   LevinIntegrator allocates the TransformationArray used by Fourier_array.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_ktGrid
use aTMDe_IO
use aTMDe_ktGrid
use aTMDe_Levin
use aTMDe_numerics
use aTMDe_control
use uTMDPDF
implicit none

!! ---- configuration ----
character(*), parameter :: inifile   = 'test.atmde'
character(*), parameter :: prefix    = 'Tests/aTMDe_ktGrid/'
real(dp),     parameter :: threshold = 1.e-4_dp

!! ---- test grid dimensions ----
integer,  parameter :: Nx   = 10
integer,  parameter :: NkT  = 40
integer,  parameter :: NQ   = 5
integer,  parameter :: Nrep = 3    !! timing repetitions

!! ---- test point ranges ----
real(dp), parameter :: xMIN_test = 0.0001_dp
!! kMIN_val must match kRanges(0) in the config *F / *p7 line
real(dp), parameter :: kMIN_val  = 0.01_dp

real(dp), parameter :: Qvalues(NQ) = [2._dp, 5._dp, 10._dp, 50._dp, 100._dp]

!! ---- arrays of test points ----
real(dp) :: xtest(Nx)
real(dp) :: kttest(NkT, NQ)

!! ---- working variables ----
integer  :: fu, ix, ikt, iq, fv, n_bad, n_bad_Q, irep
real(dp) :: time1, time2, time_make, time_extract, time_direct
real(dp), dimension(-5:5) :: from_grid, from_exact
real(dp) :: ae, re, nd, nd_min, nd_min_Q
character(len=200) :: buf

type(ktGrid)          :: grid
type(LevinIntegrator) :: levin_obj

!! shared state for toGrid_test / toFourier_inner (accessed via host association)
real(dp) :: x_cur, Q_cur
integer  :: h_cur

!! pre-computed exact values and per-(x,Q) scale for error normalisation
real(dp) :: exact_store(Nx, NkT, -5:5)
real(dp) :: f_scale_Q(Nx, -5:5)

!! ============================================================
!! Build test point arrays
!! ============================================================
do ix = 1, Nx
    xtest(ix) = 10._dp**(log10(xMIN_test) * real(Nx+1-ix, dp) / real(Nx+1, dp))
end do

do iq = 1, NQ
    do ikt = 1, NkT
        !! log scale from kMIN_val to Q; ikt=1 -> kMIN_val, ikt=NkT -> Q
        kttest(ikt, iq) = kMIN_val * (Qvalues(iq) / kMIN_val)**(real(ikt-1, dp) / real(NkT-1, dp))
    end do
end do

open(newunit=fu, file=trim(prefix)//'test.out', status='replace', action='write')

!! ============================================================
!! Initialization
!! ============================================================
call artemide_Initialize(inifile, prefix=prefix)
levin_obj = LevinIntegrator(trim(prefix)//trim(inifile), '*4   ', '*F   ', 'test_ktGrid', 1, 0)
call artemide_SetNPparameters_uTMDPDF([&
    0.874245_dp, 0.913883_dp, 0.991563_dp, 6.05412_dp, &
    0.353908_dp, 46.6064_dp,  0.115161_dp, 1.53235_dp, &
    1.31966_dp,  0.434833_dp, 0.0_dp,      0.0_dp])

!! ============================================================
!! Build ktGrid and time its construction
!! ============================================================
call cpu_time(time1)
grid = ktGrid(trim(prefix)//trim(inifile), '*4   ', '*F   ', 1, .false., 'test_ktGrid', 1)
call grid%MakeGrid(toGrid_test)
call cpu_time(time2)
time_make = time2 - time1

!! ============================================================
!! Header
!! ============================================================
call PutLine(repeat('=', 99))
call PutLine('  ktGrid accuracy test:  grid extraction vs. direct Fourier_Levin (TMD_ev_inKT)')
call PutLine('  Grid filled with: Fourier_Levin_array( uTMDPDF_inB(x,b,Q,Q^2,h) )')
call PutLine('  Reference:        uTMDPDF_inKT(x,kT,Q,Q^2,h)  [dispatches to TMD_ev_inKT]')
call PutLine('  Error metric: ae / max_kT|f(x,fv)|   (normalised by function scale per x,Q)')
call PutLine('  Flagged when normalised error > '//trim(numToStr(threshold)))
call PutLine(repeat('=', 99))

n_bad  = 0
nd_min = 1.e30_dp

!! ============================================================
!! Accuracy test — loop over Q values
!! ============================================================
do iq = 1, NQ

    call PutLine('')
    write(buf, '("  --- Q = ",F7.2," GeV  (kT from ",ES8.2," to ",F7.2,", ",I2," log-spaced points)")') &
        Qvalues(iq), kMIN_val, Qvalues(iq), NkT
    call PutLine(trim(buf))
    call PutLine('  x            kT         fl   from_grid             from_exact' &
                 //'            |error|    acc.dig')
    call PutLine('  '//repeat('-', 102))

    !! Pre-compute exact values and derive function scale for this Q
    f_scale_Q = 0._dp
    do ix = 1, Nx
        do ikt = 1, NkT
            from_exact = uTMDPDF_inKT(xtest(ix), kttest(ikt,iq), Qvalues(iq), Qvalues(iq)**2, 1)
            exact_store(ix, ikt, -5:5) = from_exact
            do fv = -5, 5
                if (fv == 0) cycle
                f_scale_Q(ix, fv) = max(f_scale_Q(ix, fv), abs(from_exact(fv)))
            end do
        end do
    end do

    n_bad_Q  = 0
    nd_min_Q = 1.e30_dp

    do ix = 1, Nx
        do ikt = 1, NkT
            from_grid = grid%Extract(xtest(ix), kttest(ikt,iq), Qvalues(iq), 1)
            do fv = -5, 5
                if (fv == 0) cycle
                ae = abs(from_grid(fv) - exact_store(ix, ikt, fv))
                re = ae / max(f_scale_Q(ix, fv), 1.e-32_dp)
                nd = -log10(max(re, 1.e-16_dp))
                if (nd < nd_min_Q) nd_min_Q = nd
                if (re > threshold) then
                    n_bad_Q = n_bad_Q + 1
                    n_bad   = n_bad + 1
                    write(buf, '(2X,ES10.3,2X,ES10.3,2X,I3,2X,ES22.14,2X,ES22.14,2X,ES10.3,2X,F7.2)') &
                        xtest(ix), kttest(ikt,iq), fv, from_grid(fv), exact_store(ix,ikt,fv), ae, nd
                    call PutLine(trim(buf))
                end if
            end do
        end do
    end do

    if (nd_min_Q < nd_min) nd_min = nd_min_Q

    write(buf, '("  Flagged: ",I5,"    Worst acc.digits: ",F7.2)') n_bad_Q, nd_min_Q
    call PutLine(trim(buf))

end do

!! ============================================================
!! Timing: grid extraction vs. direct TMD_ev_inKT  (at Q=10 GeV)
!! ============================================================
iq = 3   !! Q = 10 GeV as representative mid-range value

call cpu_time(time1)
do irep = 1, Nrep
    do ix = 1, Nx
        do ikt = 1, NkT
            from_grid = grid%Extract(xtest(ix), kttest(ikt,iq), Qvalues(iq), 1)
        end do
    end do
end do
call cpu_time(time2)
time_extract = time2 - time1

call cpu_time(time1)
do irep = 1, Nrep
    do ix = 1, Nx
        do ikt = 1, NkT
            from_exact = uTMDPDF_inKT(xtest(ix), kttest(ikt,iq), Qvalues(iq), Qvalues(iq)**2, 1)
        end do
    end do
end do
call cpu_time(time2)
time_direct = time2 - time1

!! ============================================================
!! Summary
!! ============================================================
call PutLine('')
call PutLine(repeat('=', 99))
call PutLine('  Summary')
call PutLine(repeat('=', 99))
write(buf, '("  Test points (x x kT x Q):        ",I3," x ",I3," x ",I2," = ",I6)') Nx, NkT, NQ, Nx*NkT*NQ
call PutLine(trim(buf))
write(buf, '("  Flavors tested:                  -5..5  (10 per point, total ",I7,")")') Nx*NkT*NQ*10
call PutLine(trim(buf))
write(buf, '("  Flagged cases (norm. err > ",ES8.1,"):  ",I5)') threshold, n_bad
call PutLine(trim(buf))
write(buf, '("  Worst accurate digits:           ",F7.2," (= min over all points and Q values)")') nd_min
call PutLine(trim(buf))
call PutLine('')
write(buf, '("  MakeGrid time:              ",F10.4," s")') time_make
call PutLine(trim(buf))
write(buf, '("  Grid extraction (Q=10 GeV): ",F10.4," s  (",I5," calls)")') &
    time_extract, Nrep*Nx*NkT
call PutLine(trim(buf))
write(buf, '("  Direct TMD_ev  (Q=10 GeV): ",F10.4," s  (",I5," calls)")') &
    time_direct, Nrep*Nx*NkT
call PutLine(trim(buf))
if (time_extract > 0._dp) then
    write(buf, '("  Speedup factor:             ",F10.1,"x")') time_direct / time_extract
    call PutLine(trim(buf))
end if

close(fu)

contains

!! Fill ktGrid with Fourier_Levin_array of the uTMDPDF b-space TMD at (x_in, Q_in).
!! Matches the TMDgrid_inKT interface required by ktGrid%MakeGrid.
function toGrid_test(x_in, Q_in, h_in, s1, s2)
    integer,  intent(in) :: s1, s2
    real(dp), intent(in) :: x_in, Q_in
    integer,  intent(in) :: h_in
    real(dp), dimension(1:s1, 0:s2, -5:5) :: toGrid_test
    x_cur = x_in
    Q_cur = Q_in
    h_cur = h_in
    toGrid_test = levin_obj%Fourier_array(toFourier_inner)
end function toGrid_test

!! b-space TMD at fixed (x_cur, Q_cur, Q_cur^2, h_cur).
!! This is the integrand of TMD_ev_inKT, so grid and reference share the same function.
function toFourier_inner(b)
    real(dp), intent(in) :: b
    real(dp), dimension(-5:5) :: toFourier_inner
    toFourier_inner = uTMDPDF_inB(x_cur, b, Q_cur, Q_cur**2, h_cur)
end function toFourier_inner

subroutine PutLine(str)
    character(len=*), intent(in) :: str
    write(*,  '(A)') str
    write(fu, '(A)') str
end subroutine PutLine

end program test_ktGrid
