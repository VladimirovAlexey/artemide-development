!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   test3_TMDX_DY.f90
!
!   LP reference run: computes DY cross-section at Leading Power (LP) over a
!   kinematic grid and saves results to file for later comparison with KPC.
!
!   Process: pp -> Z/gamma*  [1,1,1,3]  sqrt(s)=13 TeV  (no lepton cuts)
!     y  points: 0.0, 2.5
!     pT points: 1.5, 4.5, 11.0 GeV
!     Q  bins : [50,70], [70,90], ..., [230,250] GeV  (10 bins of 20 GeV)
!     Total   : 60 data points
!
!   NP parameters: uTMDPDF SV19 replica 0 (replace as needed)
!
!   INI file test3_LP.atmde must have:
!     *9  *p2 = F    (LP factorization)
!
!   Compile and run from the artemide root directory:
!     make program TARGET=Tests/TMDX_DY/test3_TMDX_DY.f90
!
!   Run this program FIRST; then run test3b_TMDX_DY to compare with KPC.
!   Results are saved to Tests/TMDX_DY/test3_LP.dat.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test3_TMDX_DY
use aTMDe_control
use TMDX_DY
use aTMDe_Numerics
implicit none

character(*), parameter :: inifile = 'test.atmde'
character(*), parameter :: prefix  = 'Tests/TMDX_DY/'
character(*), parameter :: datfile = 'Tests/TMDX_DY/test3_LP.dat'

!---- NP parameters (uTMDPDF SV19 replica 0) ----
real(dp), parameter :: NP(12) = [ &
    0.874245_dp, 0.913883_dp, 0.991563_dp, 6.05412_dp,  &
    0.353908_dp, 46.6064_dp,  0.115161_dp, 1.53235_dp,  &
    1.31966_dp,  0.434833_dp, 0.0_dp,      0.0_dp       ]

!---- Kinematics ----
real(dp), parameter :: s_val   = 13000._dp**2
integer,  parameter :: PROC(4) = [1, 1, 1, 3]

integer,  parameter :: NY  = 2
real(dp), parameter :: Y_LO(NY)  = [ 0._dp, 2.5_dp]
real(dp), parameter :: Y_HI(NY)  = [ 0._dp, 2.5_dp]

integer,  parameter :: NPT = 3
real(dp), parameter :: PT_LO(NPT) = [ 1.5_dp,  4.5_dp, 11._dp]
real(dp), parameter :: PT_HI(NPT) = [ 1.5_dp,  4.5_dp, 11._dp]

integer,  parameter :: NQ     = 10
real(dp), parameter :: Q_START = 50._dp
real(dp), parameter :: DQ      = 20._dp

integer,  parameter :: NTOT = NY * NPT * NQ   ! 60

integer  :: fu, iy, ipt, iq, idx
integer  :: t_cnt0, t_cnt1, t_rate
real(dp) :: xsec, Q_lo, Q_hi, pT_arr(2), Q_arr(2), y_arr(2)

!---- Initialize ----
call artemide_Initialize(inifile, prefix=prefix)
call artemide_SetNPparameters_uTMDPDF(NP)

if(.not. TMDX_DY_IsInitialized()) then
    write(*,*) 'ERROR: TMDX_DY not initialized. Check INI file (need *9 *p1 = T).'
    error stop
end if

!---- Print run header ----
write(*,'(A)') repeat('=',90)
write(*,'(A)') '  test3_TMDX_DY: LP reference  process=[1,1,1,3]  sqrt(s)=13 TeV'
write(*,'(A)') '  NP: uTMDPDF SV19 replica 0'
write(*,'(A)') repeat('=',90)
write(*,'(A)') ''
write(*,'(A)') '   iy  ipt  iq     y_lo   y_hi   pT_lo  pT_hi    Q_lo    Q_hi       xSec[pb]'
write(*,'(A)') '  '//repeat('-',80)

!---- Open output file ----
open(newunit=fu, file=datfile, status='replace', action='write')
write(fu,'(A)') '# test3_TMDX_DY: LP reference cross-sections'
write(fu,'(A)') '# process=[1,1,1,3]  sqrt(s)=13 TeV  no cuts'
write(fu,'(A)') '# NP: uTMDPDF SV19 replica 0'
write(fu,'(A)') '#'
write(fu,'(A)') '# iy  ipt  iq     y_lo   y_hi   pT_lo  pT_hi    Q_lo    Q_hi       xSec[pb]'

!---- Compute and save ----
call system_clock(t_cnt0, t_rate)
idx = 0
do iy = 1, NY
    do ipt = 1, NPT
        do iq = 1, NQ
            idx     = idx + 1
            Q_lo    = Q_START + (iq-1)*DQ
            Q_hi    = Q_START +  iq   *DQ
            pT_arr  = [PT_LO(ipt), PT_HI(ipt)]
            Q_arr   = [Q_lo, Q_hi]
            y_arr   = [Y_LO(iy), Y_HI(iy)]

            call xSec_DY(xsec, PROC, s_val, pT_arr, Q_arr, y_arr, .false.)

            write(*,  '(3I4, 2F7.3, 2F7.3, 2F8.2, ES20.12)') &
                iy, ipt, iq, Y_LO(iy), Y_HI(iy), PT_LO(ipt), PT_HI(ipt), Q_lo, Q_hi, xsec
            write(fu, '(3I4, 2F7.3, 2F7.3, 2F8.2, ES20.12)') &
                iy, ipt, iq, Y_LO(iy), Y_HI(iy), PT_LO(ipt), PT_HI(ipt), Q_lo, Q_hi, xsec
        end do
    end do
end do

call system_clock(t_cnt1)
close(fu)

write(*,'(A)') ''
write(*,'(A,I0,A)') '  Written ', NTOT, ' LP values to '//datfile
write(*,'(A,F8.3,A)') '  Evaluation time: ', &
    real(t_cnt1-t_cnt0,dp)/real(max(t_rate,1),dp), ' s'

end program test3_TMDX_DY
