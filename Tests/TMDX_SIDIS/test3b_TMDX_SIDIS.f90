!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   test3b_TMDX_SIDIS.f90
!
!   KPC comparison run: computes SIDIS cross-section with KPC factorization over the
!   same Q grid as test3_TMDX_SIDIS and compares against the saved LP results.
!
!   Process: h->pi  [1,1,1,2001]  sqrt(s)=300 GeV  (no cuts)
!     pPerp point : 0.5 GeV
!     z     point : 0.3
!     x     point : 0.1
!     Q     points: 4, 5, ..., 80 GeV  (77 points, step 1 GeV)
!
!   NP parameters: same as LP run for structural comparison;
!                  replace with KPC-fit values for physics comparison
!
!   INI file test_KPC.atmde must have:
!     *10  *p2 = T    (KPC factorization)
!
!   Compile and run from the artemide root directory:
!     make program TARGET=Tests/TMDX_SIDIS/test3b_TMDX_SIDIS.f90
!
!   Run test3_TMDX_SIDIS FIRST to generate Tests/TMDX_SIDIS/test3_LP.dat.
!   Results are written to Tests/TMDX_SIDIS/test3b_TMDX_SIDIS.out.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test3b_TMDX_SIDIS
use aTMDe_control
use TMDX_SIDIS
use aTMDe_Numerics
implicit none

character(*), parameter :: inifile = 'test_KPC.atmde'
character(*), parameter :: prefix  = 'Tests/TMDX_SIDIS/'
character(*), parameter :: datfile = 'Tests/TMDX_SIDIS/test3_LP.dat'
character(*), parameter :: outfile = 'Tests/TMDX_SIDIS/test3b_TMDX_SIDIS.out'

!---- NP parameters (replace with KPC-fit values) ----
real(dp), parameter :: NP_param(28) = [ &
    1.5_dp,       0.0859102_dp, 0.030294_dp,  0.0_dp,      &
    0.48622_dp,   0.0411753_dp, 0.569024_dp,  0.146933_dp, &
    5.26034_dp,   21.1222_dp,   7.71185_dp,   0.1565_dp,   &
    0.240061_dp,  0.0691505_dp, 1.0_dp,       1.0_dp,      &
    0.696096_dp,  0.626588_dp,  0.00331303_dp,-0.466377_dp, &
    0.88367_dp,   0.882092_dp,  1.74168_dp,   1.15036_dp,  &
    0.610318_dp, -0.101387_dp,  0.0_dp,       0.0_dp       ]

!---- Kinematics (must match test3_TMDX_SIDIS) ----
real(dp), parameter :: s_val   = 300._dp**2
integer,  parameter :: PROC(4) = [1, 1, 1, 2001]
real(dp), parameter :: CUTS(4) = [0._dp, 0._dp, 0._dp, 0._dp]

real(dp), parameter :: PPERP = 0.5_dp
real(dp), parameter :: Z_PT  = 0.3_dp
real(dp), parameter :: X_PT  = 0.1_dp

integer,  parameter :: NQ     = 77
real(dp), parameter :: Q_START = 4._dp
real(dp), parameter :: DQ      = 1._dp

real(dp) :: xLP(NQ), xKPC(NQ)
integer  :: fu_dat, fu_out
integer  :: iq, nread
real(dp) :: Q_val
real(dp) :: pPerp_arr(2), z_arr(2), x_arr(2), Q_arr(2)
real(dp) :: rd, rd_max, t_kpc
integer  :: t_cnt0, t_cnt1, t_rate
integer  :: iq_r
character(len=200) :: line

!---- Initialize ----
call artemide_Initialize(inifile, prefix=prefix)
call artemide_SetNPparameters(NP_param)

if(.not. TMDX_SIDIS_IsInitialized()) then
    write(*,*) 'ERROR: TMDX_SIDIS not initialized. Check INI file (need *10 *p1 = T).'
    error stop
end if

!---- Fixed kinematic arrays (lo=hi → point evaluation) ----
pPerp_arr = [PPERP, PPERP]
z_arr     = [Z_PT,  Z_PT ]
x_arr     = [X_PT,  X_PT ]

!---- Read LP reference from file ----
open(newunit=fu_dat, file=datfile, status='old', action='read', err=901)
nread = 0
do
    read(fu_dat, '(A)', end=10) line
    if(line(1:1) == '#') cycle
    nread = nread + 1
    if(nread > NQ) exit
    read(line, *) iq_r, Q_val, xLP(nread)
end do
10 continue
close(fu_dat)

if(nread /= NQ) then
    write(*,'(A,I0,A,I0)') 'ERROR: expected ', NQ, ' data lines in LP file, got ', nread
    error stop
end if

!---- Compute KPC ----
call system_clock(t_cnt0, t_rate)
do iq = 1, NQ
    Q_val = Q_START + (iq-1)*DQ
    Q_arr = [Q_val, Q_val]
    call xSec_SIDIS(xKPC(iq), PROC, s_val, pPerp_arr, z_arr, x_arr, Q_arr, .false., CUTS)
end do
call system_clock(t_cnt1)
t_kpc = real(t_cnt1-t_cnt0,dp)/real(max(t_rate,1),dp)

!---- Print comparison ----
open(newunit=fu_out, file=outfile, status='replace', action='write')
call PrintHeader()

rd_max = 0._dp
do iq = 1, NQ
    Q_val = Q_START + (iq-1)*DQ
    rd    = abs(xKPC(iq) - xLP(iq)) / max(abs(xLP(iq)), 1.e-300_dp)
    rd_max = max(rd_max, rd)
    call PrintRow(iq, Q_val, xLP(iq), xKPC(iq), rd)
end do

call PrintSummary(rd_max, t_kpc)

close(fu_out)
stop

901 write(*,'(A)') 'ERROR: cannot open LP reference file: '//datfile
    write(*,'(A)') '       Run test3_TMDX_SIDIS first.'
    error stop

!=======================================================================
contains
!=======================================================================

subroutine OutLine(str)
    character(len=*), intent(in) :: str
    write(*,'(A)') str
    write(fu_out,'(A)') str
end subroutine OutLine

subroutine PrintHeader()
    call OutLine(repeat('=',100))
    call OutLine('  test3b_TMDX_SIDIS: KPC vs LP comparison')
    call OutLine('  process=[1,1,1,2001]  sqrt(s)=300 GeV  no cuts')
    write(*,'(A,F4.1,A,F4.1,A,F4.1)') '  pPerp=', PPERP, ' GeV   z=', Z_PT, '   x=', X_PT
    write(fu_out,'(A,F4.1,A,F4.1,A,F4.1)') '  pPerp=', PPERP, ' GeV   z=', Z_PT, '   x=', X_PT
    call OutLine('  NP: same for LP and KPC (replace with KPC-fit values for physics)')
    call OutLine('  LP  factorization: test.atmde      (*10 *p2 = F)')
    call OutLine('  KPC factorization: test_KPC.atmde  (*10 *p2 = T)')
    call OutLine('  rd = |KPC - LP| / |LP|')
    call OutLine(repeat('=',100))
    call OutLine('')
    call OutLine('     iq     Q [GeV]         LP [pb]             KPC [pb]              rd')
    call OutLine('  '//repeat('-',80))
end subroutine PrintHeader

subroutine PrintRow(iq, Q_val, lp, kpc, rd)
    integer,  intent(in) :: iq
    real(dp), intent(in) :: Q_val, lp, kpc, rd
    character(len=200) :: ln
    write(ln,'(I6, F12.3, 2X, ES18.9, 3X, ES18.9, 3X, ES10.3)') iq, Q_val, lp, kpc, rd
    call OutLine(trim(ln))
end subroutine PrintRow

subroutine PrintSummary(rd_max, t_kpc)
    real(dp), intent(in) :: rd_max, t_kpc
    character(len=100) :: ln
    call OutLine('')
    call OutLine(repeat('=',100))
    write(ln,'(2X,"Overall max rd (KPC vs LP) across all ",I0," Q points: ",ES10.3)') NQ, rd_max
    call OutLine(trim(ln))
    write(ln,'(2X,"Evaluation time (KPC): ",F8.3," s")') t_kpc
    call OutLine(trim(ln))
    call OutLine(repeat('=',100))
end subroutine PrintSummary

end program test3b_TMDX_SIDIS
