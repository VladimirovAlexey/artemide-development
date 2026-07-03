!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   test3_TMDX_SIDIS.f90
!
!   LP reference run: computes SIDIS cross-section at Leading Power (LP) over a
!   Q grid and saves results to file for later comparison with KPC.
!
!   Process: h->pi  [1,1,1,2001]  sqrt(s)=300 GeV  (no cuts)
!     pPerp point : 0.5 GeV
!     z     point : 0.3
!     x     point : 0.1
!     Q     points: 4, 5, ..., 80 GeV  (77 points, step 1 GeV)
!
!   NP parameters: replace with fit values as needed
!
!   INI file test.atmde must have:
!     *10  *p2 = F    (LP factorization)
!
!   Compile and run from the artemide root directory:
!     make program TARGET=Tests/TMDX_SIDIS/test3_TMDX_SIDIS.f90
!
!   Run this program FIRST; then run test3b_TMDX_SIDIS to compare with KPC.
!   Results are saved to Tests/TMDX_SIDIS/test3_LP.dat.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test3_TMDX_SIDIS
use aTMDe_control
use TMDX_SIDIS
use aTMDe_Numerics
implicit none

character(*), parameter :: inifile = 'test.atmde'
character(*), parameter :: prefix  = 'Tests/TMDX_SIDIS/'
character(*), parameter :: datfile = 'Tests/TMDX_SIDIS/test3_LP.dat'

!---- NP parameters (replace with fit values) ----
real(dp), parameter :: NP_param(28) = [ &
    1.5_dp,       0.0859102_dp, 0.030294_dp,  0.0_dp,      &
    0.48622_dp,   0.0411753_dp, 0.569024_dp,  0.146933_dp, &
    5.26034_dp,   21.1222_dp,   7.71185_dp,   0.1565_dp,   &
    0.240061_dp,  0.0691505_dp, 1.0_dp,       1.0_dp,      &
    0.696096_dp,  0.626588_dp,  0.00331303_dp,-0.466377_dp, &
    0.88367_dp,   0.882092_dp,  1.74168_dp,   1.15036_dp,  &
    0.610318_dp, -0.101387_dp,  0.0_dp,       0.0_dp       ]

!---- Kinematics ----
real(dp), parameter :: s_val   = 300._dp**2
integer,  parameter :: PROC(4) = [1, 1, 1, 2001]
real(dp), parameter :: CUTS(4) = [0._dp, 0._dp, 0._dp, 0._dp]

real(dp), parameter :: PPERP = 0.5_dp
real(dp), parameter :: Z_PT  = 0.3_dp
real(dp), parameter :: X_PT  = 0.1_dp

integer,  parameter :: NQ     = 77
real(dp), parameter :: Q_START = 4._dp
real(dp), parameter :: DQ      = 1._dp

integer  :: fu, iq
integer  :: t_cnt0, t_cnt1, t_rate
real(dp) :: xsec, Q_val
real(dp) :: pPerp_arr(2), z_arr(2), x_arr(2), Q_arr(2)

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

!---- Print run header ----
write(*,'(A)') repeat('=',90)
write(*,'(A)') '  test3_TMDX_SIDIS: LP reference  process=[1,1,1,2001]  sqrt(s)=300 GeV'
write(*,'(A,F4.1,A,F4.1,A,F4.1)') '  pPerp=', PPERP, ' GeV   z=', Z_PT, '   x=', X_PT
write(*,'(A)') repeat('=',90)
write(*,'(A)') ''
write(*,'(A)') '     iq     Q [GeV]       xSec [pb]'
write(*,'(A)') '  '//repeat('-',45)

!---- Open output file ----
open(newunit=fu, file=datfile, status='replace', action='write')
write(fu,'(A)') '# test3_TMDX_SIDIS: LP reference cross-sections'
write(fu,'(A)') '# process=[1,1,1,2001]  sqrt(s)=300 GeV  no cuts'
write(fu,'(A,F4.1,A,F4.1,A,F4.1)') '# pPerp=', PPERP, ' GeV   z=', Z_PT, '   x=', X_PT
write(fu,'(A)') '#'
write(fu,'(A)') '#  iq     Q [GeV]       xSec [pb]'

!---- Compute and save ----
call system_clock(t_cnt0, t_rate)
do iq = 1, NQ
    Q_val = Q_START + (iq-1)*DQ
    Q_arr = [Q_val, Q_val]

    call xSec_SIDIS(xsec, PROC, s_val, pPerp_arr, z_arr, x_arr, Q_arr, .false., CUTS)

    write(*,  '(I6, F12.3, ES20.12)') iq, Q_val, xsec
    write(fu, '(I6, F12.3, ES20.12)') iq, Q_val, xsec
end do
call system_clock(t_cnt1)

close(fu)

write(*,'(A)') ''
write(*,'(A,I0,A)') '  Written ', NQ, ' LP values to '//datfile
write(*,'(A,F8.3,A)') '  Evaluation time: ', &
    real(t_cnt1-t_cnt0,dp)/real(max(t_rate,1),dp), ' s'

end program test3_TMDX_SIDIS
