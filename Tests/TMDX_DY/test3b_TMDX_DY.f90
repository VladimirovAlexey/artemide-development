!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   test3b_TMDX_DY.f90
!
!   KPC comparison run: computes DY cross-section with KPC factorization over the
!   same kinematic grid as test3_TMDX_DY and compares against the saved LP results.
!
!   Process: pp -> Z/gamma*  [1,1,1,3]  sqrt(s)=13 TeV  (no lepton cuts)
!     y  points: 0.0, 2.5
!     pT points: 1.5, 4.5, 11.0 GeV
!     Q  bins : [50,70], [70,90], ..., [230,250] GeV  (10 bins of 20 GeV)
!     Total   : 60 data points
!
!   NP parameters: uTMDPDF SV19 replica 0 (same as LP run for structural comparison;
!                  replace with KPC-fit NP values for physics comparison)
!
!   INI file test3_KPC.atmde must have:
!     *9  *p2 = T    (KPC factorization)
!
!   Compile and run from the artemide root directory:
!     make program TARGET=Tests/TMDX_DY/test3b_TMDX_DY.f90
!
!   Run test3_TMDX_DY FIRST to generate Tests/TMDX_DY/test3_LP.dat.
!   Results are written to Tests/TMDX_DY/test3b_TMDX_DY.out.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test3b_TMDX_DY
use aTMDe_control
use TMDX_DY
use aTMDe_Numerics
implicit none

character(*), parameter :: inifile = 'test_KPC.atmde'
character(*), parameter :: prefix  = 'Tests/TMDX_DY/'
character(*), parameter :: datfile = 'Tests/TMDX_DY/test3_LP.dat'
character(*), parameter :: outfile = 'Tests/TMDX_DY/test3b_TMDX_DY.out'

!---- NP parameters (uTMDPDF SV19 replica 0; replace with KPC-fit values) ----
real(dp), parameter :: NP(12) = [ &
    0.874245_dp, 0.913883_dp, 0.991563_dp, 6.05412_dp,  &
    0.353908_dp, 46.6064_dp,  0.115161_dp, 1.53235_dp,  &
    1.31966_dp,  0.434833_dp, 0.0_dp,      0.0_dp       ]

!---- Kinematics (must match test3_TMDX_DY) ----
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

real(dp) :: xLP(NTOT), xKPC(NTOT)
integer  :: fu_dat, fu_out
integer  :: iy, ipt, iq, idx, nread
integer  :: iy_r, ipt_r, iq_r
real(dp) :: y_lo_r, y_hi_r, pt_lo_r, pt_hi_r, q_lo_r, q_hi_r
real(dp) :: Q_lo, Q_hi, pT_arr(2), Q_arr(2), y_arr(2)
real(dp) :: rd, rd_max, t_kpc
integer  :: t_cnt0, t_cnt1, t_rate
character(len=300) :: line

!---- Initialize ----
call artemide_Initialize(inifile, prefix=prefix)
call artemide_SetNPparameters_uTMDPDF(NP)

if(.not. TMDX_DY_IsInitialized()) then
    write(*,*) 'ERROR: TMDX_DY not initialized. Check INI file (need *9 *p1 = T).'
    error stop
end if

!---- Read LP reference from file ----
open(newunit=fu_dat, file=datfile, status='old', action='read', err=901)
nread = 0
do
    read(fu_dat, '(A)', end=10) line
    if(line(1:1) == '#') cycle
    nread = nread + 1
    if(nread > NTOT) exit
    read(line, *) iy_r, ipt_r, iq_r, &
                  y_lo_r, y_hi_r, pt_lo_r, pt_hi_r, q_lo_r, q_hi_r, &
                  xLP(nread)
end do
10 continue
close(fu_dat)

if(nread /= NTOT) then
    write(*,'(A,I0,A,I0)') 'ERROR: expected ', NTOT, ' data lines in LP file, got ', nread
    error stop
end if

!---- Compute KPC ----
call system_clock(t_cnt0, t_rate)
idx = 0
do iy = 1, NY
    do ipt = 1, NPT
        do iq = 1, NQ
            idx    = idx + 1
            Q_lo   = Q_START + (iq-1)*DQ
            Q_hi   = Q_START +  iq   *DQ
            pT_arr = [PT_LO(ipt), PT_HI(ipt)]
            Q_arr  = [Q_lo, Q_hi]
            y_arr  = [Y_LO(iy), Y_HI(iy)]
            call xSec_DY(xKPC(idx), PROC, s_val, pT_arr, Q_arr, y_arr, .false.)
        end do
    end do
end do

call system_clock(t_cnt1)
t_kpc = real(t_cnt1-t_cnt0,dp)/real(max(t_rate,1),dp)

!---- Print comparison ----
open(newunit=fu_out, file=outfile, status='replace', action='write')
call PrintHeader()

rd_max = 0._dp
do iy = 1, NY
    do ipt = 1, NPT
        call SectionHeader(iy, ipt)
        do iq = 1, NQ
            idx = (iy-1)*NPT*NQ + (ipt-1)*NQ + iq
            rd  = abs(xKPC(idx) - xLP(idx)) / max(abs(xLP(idx)), 1.e-300_dp)
            rd_max = max(rd_max, rd)
            call PrintRow(iq, xLP(idx), xKPC(idx), rd)
        end do
    end do
end do

call PrintSummary(rd_max, t_kpc)

close(fu_out)
stop

901 write(*,'(A)') 'ERROR: cannot open LP reference file: '//datfile
    write(*,'(A)') '       Run test3_TMDX_DY first.'
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
    call OutLine('  test3b_TMDX_DY: KPC vs LP comparison')
    call OutLine('  process=[1,1,1,3]  sqrt(s)=13 TeV  no cuts')
    call OutLine('  NP: uTMDPDF SV19 replica 0  (same for LP and KPC)')
    call OutLine('  LP  factorization: test3_LP.atmde   (*9 *p2 = F)')
    call OutLine('  KPC factorization: test3_KPC.atmde  (*9 *p2 = T)')
    call OutLine('  rd = |KPC - LP| / |LP|')
    call OutLine(repeat('=',100))
end subroutine PrintHeader

subroutine SectionHeader(iy, ipt)
    integer, intent(in) :: iy, ipt
    character(len=120) :: ln
    call OutLine('')
    call OutLine(repeat('-',80))
    write(ln,'(2X,"y bin: [",F6.3,",",F6.3,"]   pT bin: [",F5.1,",",F5.1,"] GeV")') &
        Y_LO(iy), Y_HI(iy), PT_LO(ipt), PT_HI(ipt)
    call OutLine(trim(ln))
    call OutLine('  Q range [GeV]          LP [pb]             KPC [pb]              rd')
    call OutLine('  '//repeat('-',78))
end subroutine SectionHeader

subroutine PrintRow(iq, lp, kpc, rd)
    integer,  intent(in) :: iq
    real(dp), intent(in) :: lp, kpc, rd
    real(dp) :: q_lo, q_hi
    character(len=200) :: ln
    q_lo = Q_START + (iq-1)*DQ
    q_hi = Q_START +  iq   *DQ
    write(ln,'(2X,"[",F6.1,",",F6.1,"]",2X,ES18.9,3X,ES18.9,3X,ES10.3)') &
        q_lo, q_hi, lp, kpc, rd
    call OutLine(trim(ln))
end subroutine PrintRow

subroutine PrintSummary(rd_max, t_kpc)
    real(dp), intent(in) :: rd_max, t_kpc
    character(len=100) :: ln
    call OutLine('')
    call OutLine(repeat('=',100))
    write(ln,'(2X,"Overall max rd (KPC vs LP) across all ",I0," bins: ",ES10.3)') NTOT, rd_max
    call OutLine(trim(ln))
    write(ln,'(2X,"Evaluation time (KPC): ",F8.3," s")') t_kpc
    call OutLine(trim(ln))
    call OutLine(repeat('=',100))
end subroutine PrintSummary

end program test3b_TMDX_DY
