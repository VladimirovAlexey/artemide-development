!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   test_TMDX_DY.f90
!
!   Benchmark and comparison of three qT-spectrum evaluation methods:
!     (1) xSec_DY       -- single-bin call, sequential loop
!     (2) xSec_DY_List  -- list call, no qT-partitioning
!     (3) xSec_DY_List  -- list call, with qT-partitioning (Chebyshev interpolation)
!
!   Case A: Z-boson production, pp @ sqrt(s)=13 TeV
!     process = [1,1,1,3],  Q in [70,110] GeV,  no lepton cuts
!     y bins  : [-1,1], [1,2], [2,3]
!     qT range: 0--20 GeV
!       fine bins  : 40 bins of dqT=0.5 GeV
!       coarse bins: 10 bins of dqT=2.0 GeV
!
!   Case B: Drell-Yan at sqrt(s)=39 GeV (E288-like), process=[2,1,1,1]
!     xF in [-0.1, 0.5],  Q bins [11,12], [12,14], [14,16] GeV
!     qT range: 0--2 GeV,  10 bins of dqT=0.2 GeV
!
!   NP parameters: uTMDPDF SV19 replica 0
!
!   Compile and run from the artemide root directory:
!     make program TARGET=Tests/TMDX_DY/test_TMDX_DY.f90
!
!   The INI file must be present at Tests/TMDX_DY/test.atmde.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_TMDX_DY
use aTMDe_control
use TMDX_DY
use aTMDe_Numerics
implicit none

!---- paths ----
character(*), parameter :: inifile = 'test.atmde'
character(*), parameter :: prefix  = 'Tests/TMDX_DY/'
character(*), parameter :: outfile = 'Tests/TMDX_DY/test_TMDX_DY.out'

!---- NP parameters (uTMDPDF SV19 replica 0) ----
real(dp), parameter :: NP(12) = [ &
    0.874245_dp, 0.913883_dp, 0.991563_dp, 6.05412_dp,  &
    0.353908_dp, 46.6064_dp,  0.115161_dp, 1.53235_dp,  &
    1.31966_dp,  0.434833_dp, 0.0_dp,      0.0_dp       ]

!=======================================================================
! Case A: Z production, pp @ 13 TeV
!=======================================================================
integer,  parameter :: PROC_Z(4) = [1, 1, 1, 3]
real(dp), parameter :: s_Z       = 13000._dp**2
real(dp), parameter :: QMIN_Z    = 70._dp
real(dp), parameter :: QMAX_Z    = 110._dp

!---- y bins: YBINS(iy,1)=ymin, YBINS(iy,2)=ymax  (column-major reshape) ----
integer,  parameter :: NY = 3
real(dp), parameter :: YBINS(NY,2) = reshape( &
    [-1._dp, 1._dp, 2._dp,   1._dp, 2._dp, 3._dp], [NY,2])

integer,  parameter :: NQT_F = 40,  NQT_C = 10
real(dp), parameter :: DQT_F = 0.5_dp, DQT_C = 2.0_dp

integer, parameter :: NF = NY * NQT_F   ! 120 bins (fine)
integer, parameter :: NC = NY * NQT_C   !  30 bins (coarse)

real(dp) :: XF_1pt(NF), XF_list(NF), XF_part(NF)
real(dp) :: qT_F(NF,2), Q_F(NF,2), y_F(NF,2), sArr_F(NF), cutP_F(NF,4)
integer  :: proc_F(NF,4)
logical  :: cuts_F(NF)

real(dp) :: XC_1pt(NC), XC_list(NC), XC_part(NC)
real(dp) :: qT_C(NC,2), Q_C(NC,2), y_C(NC,2), sArr_C(NC), cutP_C(NC,4)
integer  :: proc_C(NC,4)
logical  :: cuts_C(NC)

real(dp) :: tF_1pt, tF_list, tF_part
real(dp) :: tC_1pt, tC_list, tC_part

!=======================================================================
! Case B: E288-like, sqrt(s)=39 GeV, process=[2,1,1,1]
!=======================================================================
integer,  parameter :: PROC_E_DEF(4) = [2, 1, 1, 1]
real(dp), parameter :: s_E288     = 39._dp**2
real(dp), parameter :: XF_E_MIN   = -0.1_dp
real(dp), parameter :: XF_E_MAX   =  0.5_dp
integer,  parameter :: NQ_E       = 3
integer,  parameter :: NQT_E      = 10
real(dp), parameter :: DQT_E      = 0.2_dp

real(dp), parameter :: QBINS_E(NQ_E,2) = reshape( &
    [11._dp, 12._dp, 14._dp,   12._dp, 14._dp, 16._dp], [NQ_E,2])

integer, parameter :: NE = NQ_E * NQT_E   ! 30 bins

real(dp) :: XE_1pt(NE), XE_list(NE), XE_part(NE)
real(dp) :: qT_E(NE,2), Q_E(NE,2), y_E(NE,2), sArr_E(NE), cutP_E(NE,4)
integer  :: proc_E(NE,4)
logical  :: cuts_E(NE)

real(dp) :: tE_1pt, tE_list, tE_part

!=======================================================================
! Accumulated precision stats from each PrintTableGen call
! (3 calls: A-fine, A-coarse, B)
!=======================================================================
real(dp) :: prec_avg_list(3), prec_max_list(3)
real(dp) :: prec_avg_part(3), prec_max_part(3)

integer :: fu, iy, iq, iqt, idx

!======================================================================
! Initialize
!======================================================================
call artemide_Initialize(inifile, prefix=prefix)
call artemide_SetNPparameters_uTMDPDF(NP)

if(.not. TMDX_DY_IsInitialized()) then
    write(*,*) 'ERROR: TMDX_DY is not initialized. Enable it in the ini file.'
    error stop
end if

open(newunit=fu, file=outfile, status='replace', action='write')

!======================================================================
! Build Case A (Z) list arrays
!======================================================================
do iy = 1, NY
    do iqt = 1, NQT_F
        idx = (iy-1)*NQT_F + iqt
        qT_F(idx,1)    = (iqt-1)*DQT_F
        qT_F(idx,2)    =  iqt   *DQT_F
        Q_F(idx,1)     = QMIN_Z;   Q_F(idx,2)  = QMAX_Z
        y_F(idx,1)     = YBINS(iy,1); y_F(idx,2) = YBINS(iy,2)
        sArr_F(idx)    = s_Z
        proc_F(idx,:)  = PROC_Z
        cuts_F(idx)    = .false.;  cutP_F(idx,:) = 0._dp
    end do
end do

do iy = 1, NY
    do iqt = 1, NQT_C
        idx = (iy-1)*NQT_C + iqt
        qT_C(idx,1)    = (iqt-1)*DQT_C
        qT_C(idx,2)    =  iqt   *DQT_C
        Q_C(idx,1)     = QMIN_Z;   Q_C(idx,2)  = QMAX_Z
        y_C(idx,1)     = YBINS(iy,1); y_C(idx,2) = YBINS(iy,2)
        sArr_C(idx)    = s_Z
        proc_C(idx,:)  = PROC_Z
        cuts_C(idx)    = .false.;  cutP_C(idx,:) = 0._dp
    end do
end do

!======================================================================
! Build Case B (E288) list arrays
! For proc1=2 the y-array carries the xF integration limits
!======================================================================
do iq = 1, NQ_E
    do iqt = 1, NQT_E
        idx = (iq-1)*NQT_E + iqt
        qT_E(idx,1)    = (iqt-1)*DQT_E
        qT_E(idx,2)    =  iqt   *DQT_E
        Q_E(idx,1)     = QBINS_E(iq,1); Q_E(idx,2) = QBINS_E(iq,2)
        y_E(idx,1)     = XF_E_MIN       ! xF limits via y-array (proc1=2)
        y_E(idx,2)     = XF_E_MAX
        sArr_E(idx)    = s_E288
        proc_E(idx,:)  = PROC_E_DEF
        cuts_E(idx)    = .false.;  cutP_E(idx,:) = 0._dp
    end do
end do

call PrintHeader()

!======================================================================
! CASE A -- FINE BINS: dqT=0.5 GeV, 40 per y-bin, 120 total
!======================================================================
call SectionHeader('CASE A -- FINE BINS: dqT=0.5 GeV, N_qT=40 per y-bin, 120 bins total')

call system_clock_start(tF_1pt)
do idx = 1, NF
    call xSec_DY(XF_1pt(idx), proc_F(idx,:), sArr_F(idx), &
                 qT_F(idx,:), Q_F(idx,:), y_F(idx,:), .false.)
end do
call system_clock_stop(tF_1pt)

call system_clock_start(tF_list)
call xSec_DY_List(XF_list, proc_F, sArr_F, qT_F, Q_F, y_F, cuts_F, cutP_F, &
                  doPartitioning=.false.)
call system_clock_stop(tF_list)

call system_clock_start(tF_part)
call xSec_DY_List(XF_part, proc_F, sArr_F, qT_F, Q_F, y_F, cuts_F, cutP_F, &
                  doPartitioning=.true.)
call system_clock_stop(tF_part)

call PrintTimings(tF_1pt, tF_list, tF_part)
call PrintTableGen(NY, NQT_F, NF, DQT_F, YBINS(:,1), YBINS(:,2), 'y bin    ', &
                   XF_1pt, XF_list, XF_part, &
                   prec_avg_list(1), prec_max_list(1), prec_avg_part(1), prec_max_part(1))

!======================================================================
! CASE A -- COARSE BINS: dqT=2 GeV, 10 per y-bin, 30 total
!======================================================================
call SectionHeader('CASE A -- COARSE BINS: dqT=2 GeV, N_qT=10 per y-bin, 30 bins total')

call system_clock_start(tC_1pt)
do idx = 1, NC
    call xSec_DY(XC_1pt(idx), proc_C(idx,:), sArr_C(idx), &
                 qT_C(idx,:), Q_C(idx,:), y_C(idx,:), .false.)
end do
call system_clock_stop(tC_1pt)

call system_clock_start(tC_list)
call xSec_DY_List(XC_list, proc_C, sArr_C, qT_C, Q_C, y_C, cuts_C, cutP_C, &
                  doPartitioning=.false.)
call system_clock_stop(tC_list)

call system_clock_start(tC_part)
call xSec_DY_List(XC_part, proc_C, sArr_C, qT_C, Q_C, y_C, cuts_C, cutP_C, &
                  doPartitioning=.true.)
call system_clock_stop(tC_part)

call PrintTimings(tC_1pt, tC_list, tC_part)
call PrintTableGen(NY, NQT_C, NC, DQT_C, YBINS(:,1), YBINS(:,2), 'y bin    ', &
                   XC_1pt, XC_list, XC_part, &
                   prec_avg_list(2), prec_max_list(2), prec_avg_part(2), prec_max_part(2))

!======================================================================
! CASE B -- E288-like: dqT=0.2 GeV, 10 per Q-bin, 30 total
!======================================================================
call SectionHeader('CASE B -- E288-like: dqT=0.2 GeV, N_qT=10 per Q-bin, 30 bins total')
call OutLine('  process=[2,1,1,1]  sqrt(s)=39 GeV  xF in [-0.1, 0.5]')
call OutLine('  Q bins: [11,12], [12,14], [14,16] GeV   qT: 0--2 GeV')

call system_clock_start(tE_1pt)
do idx = 1, NE
    call xSec_DY(XE_1pt(idx), proc_E(idx,:), sArr_E(idx), &
                 qT_E(idx,:), Q_E(idx,:), y_E(idx,:), .false.)
end do
call system_clock_stop(tE_1pt)

call system_clock_start(tE_list)
call xSec_DY_List(XE_list, proc_E, sArr_E, qT_E, Q_E, y_E, cuts_E, cutP_E, &
                  doPartitioning=.false.)
call system_clock_stop(tE_list)

call system_clock_start(tE_part)
call xSec_DY_List(XE_part, proc_E, sArr_E, qT_E, Q_E, y_E, cuts_E, cutP_E, &
                  doPartitioning=.true.)
call system_clock_stop(tE_part)

call PrintTimings(tE_1pt, tE_list, tE_part)
call PrintTableGen(NQ_E, NQT_E, NE, DQT_E, QBINS_E(:,1), QBINS_E(:,2), 'Q bin    ', &
                   XE_1pt, XE_list, XE_part, &
                   prec_avg_list(3), prec_max_list(3), prec_avg_part(3), prec_max_part(3))

call PrintTotalSummary()
call PrintTotalPrecisionSummary()

close(fu)

!=======================================================================
contains
!=======================================================================

subroutine system_clock_start(t)
    real(dp), intent(out) :: t
    integer :: cnt, rate
    call system_clock(cnt, rate)
    t = real(cnt, dp) / real(max(rate,1), dp)
end subroutine system_clock_start

subroutine system_clock_stop(t)
    real(dp), intent(inout) :: t
    integer :: cnt, rate
    call system_clock(cnt, rate)
    t = real(cnt, dp) / real(max(rate,1), dp) - t
end subroutine system_clock_stop

!-----------------------------------------------------------------------

subroutine OutLine(str)
    character(len=*), intent(in) :: str
    write(*,'(A)') str
    write(fu,'(A)') str
end subroutine OutLine

subroutine SectionHeader(title)
    character(len=*), intent(in) :: title
    call OutLine('')
    call OutLine(repeat('=',100))
    call OutLine('  '//trim(title))
    call OutLine(repeat('=',100))
end subroutine SectionHeader

subroutine PrintHeader()
    call OutLine(repeat('=',100))
    call OutLine('  TMDX_DY test')
    call OutLine('  NP (uTMDPDF SV19 rep.0):  0.874245  0.913883  0.991563  6.05412')
    call OutLine('                              0.353908  46.6064   0.115161  1.53235')
    call OutLine('                              1.31966   0.434833  0         0')
    call OutLine('  Case A: Z production, pp @ 13 TeV, process=[1,1,1,3], Q in [70,110] GeV')
    call OutLine('           y bins: [-1,1],[1,2],[2,3]   qT: 0--20 GeV, no lepton cuts')
    call OutLine('  Case B: DY at sqrt(s)=39 GeV, process=[2,1,1,1], xF in [-0.1,0.5]')
    call OutLine('           Q bins: [11,12],[12,14],[14,16] GeV   qT: 0--2 GeV, dqT=0.2 GeV')
    call OutLine(repeat('=',100))
end subroutine PrintHeader

subroutine PrintTimings(t1pt, tlist, tpart)
    real(dp), intent(in) :: t1pt, tlist, tpart
    character(len=80) :: ln
    call OutLine('')
    call OutLine('  Timings (wall clock):')
    write(ln,'(4X,"single-point loop  : ",F8.3," s")') t1pt;   call OutLine(trim(ln))
    write(ln,'(4X,"list (no partition): ",F8.3," s")') tlist;  call OutLine(trim(ln))
    write(ln,'(4X,"list (partitioned) : ",F8.3," s")') tpart;  call OutLine(trim(ln))
    write(ln,'(4X,"speedup list / 1pt : ",F7.2)') t1pt / max(tlist, 1.e-9_dp)
    call OutLine(trim(ln))
    write(ln,'(4X,"speedup part / 1pt : ",F7.2)') t1pt / max(tpart, 1.e-9_dp)
    call OutLine(trim(ln))
end subroutine PrintTimings

!  General table printer.
!
!  Column layout (1-indexed character positions in each output line):
!    1-2   : 2X indent
!    3-17  : qT bin "[F5.2,F5.2]  " or stats label (A15)
!    18-32 : group bin "[F6.2,F6.2]"
!    33-35 : 3X gap
!    36-49 : 1-point value (ES14.6) or 14 dashes
!    50-52 : 3X gap
!    53-66 : list value (ES14.6) or 14 dashes
!    67-69 : 3X gap
!    70-78 : rd(list) (ES9.2) -- avg or max for stats rows
!    79-81 : 3X gap
!    82-95 : partitioned value (ES14.6) or 14 dashes
!    96-98 : 3X gap
!    99-107: rd(part) (ES9.2) -- avg or max for stats rows
!
!  After each group's data rows, two stats lines are printed (avg rd, max rd)
!  with 14 dashes in all cross-section value columns.
!
!  rd = |X_method - X_1pt| / |X_1pt|  (X_1pt = single-point reference)
!
!  Optional out arguments return overall avg and max rd across all bins.
subroutine PrintTableGen(ngroup, nqt, ntot, dqt, grp_lo, grp_hi, grp_hdr, X1, X2, X3, &
                         ret_avg_rd_list, ret_max_rd_list, ret_avg_rd_part, ret_max_rd_part)
    integer,          intent(in)  :: ngroup, nqt, ntot
    real(dp),         intent(in)  :: dqt
    real(dp),         intent(in)  :: grp_lo(ngroup), grp_hi(ngroup)
    character(len=*), intent(in)  :: grp_hdr
    real(dp),         intent(in)  :: X1(ntot), X2(ntot), X3(ntot)
    real(dp), intent(out), optional :: ret_avg_rd_list, ret_max_rd_list
    real(dp), intent(out), optional :: ret_avg_rd_part, ret_max_rd_part

    integer  :: ig, iqt, idx
    real(dp) :: rd2, rd3
    real(dp) :: sum2_grp, max2_grp, sum3_grp, max3_grp
    real(dp) :: sum2_all,  max2_all,  sum3_all,  max3_all
    character(len=200) :: ln

    ! 15-char label for stats rows; 14-char dash string for X-value columns
    character(len=15), parameter :: LBL_AVG = '[avg rd]       '
    character(len=15), parameter :: LBL_MAX = '[max rd]       '
    character(len=14), parameter :: DASH14  = '--------------'

    ! format strings
    character(*), parameter :: FMT_DATA = &
      '(2X,"[",F5.2,",",F5.2,"]  ","[",F6.2,",",F6.2,"]",'// &
      '3X,ES14.6,3X,ES14.6,3X,ES9.2,3X,ES14.6,3X,ES9.2)'
    character(*), parameter :: FMT_STAT = &
      '(2X,A15,"[",F6.2,",",F6.2,"]",'// &
      '3X,A14,3X,A14,3X,ES9.2,3X,A14,3X,ES9.2)'

    sum2_all = 0._dp;  max2_all = 0._dp
    sum3_all = 0._dp;  max3_all = 0._dp

    call OutLine('')
    call OutLine( &
      '  rd = |X_method - X_1pt| / |X_1pt|  (X_1pt = single-point sequential reference)')
    call OutLine( &
      '  qT bin [GeV]    '//trim(grp_hdr)//'       '// &
      '  1-point [pb]      list [pb]    rd(list)'// &
      '    partitioned [pb]   rd(part)')
    call OutLine('  '//repeat('-',125))

    do ig = 1, ngroup
        sum2_grp = 0._dp;  max2_grp = 0._dp
        sum3_grp = 0._dp;  max3_grp = 0._dp

        do iqt = 1, nqt
            idx = (ig-1)*nqt + iqt
            rd2 = abs(X2(idx)-X1(idx)) / max(abs(X1(idx)), 1.e-300_dp)
            rd3 = abs(X3(idx)-X1(idx)) / max(abs(X1(idx)), 1.e-300_dp)
            write(ln, FMT_DATA) &
                (iqt-1)*dqt, iqt*dqt, grp_lo(ig), grp_hi(ig), &
                X1(idx), X2(idx), rd2, X3(idx), rd3
            call OutLine(trim(ln))
            sum2_grp = sum2_grp + rd2;  max2_grp = max(max2_grp, rd2)
            sum3_grp = sum3_grp + rd3;  max3_grp = max(max3_grp, rd3)
        end do

        ! stats lines: dashes in X columns, avg/max rd in rd columns
        write(ln, FMT_STAT) LBL_AVG, grp_lo(ig), grp_hi(ig), &
            DASH14, DASH14, sum2_grp/nqt, DASH14, sum3_grp/nqt
        call OutLine(trim(ln))
        write(ln, FMT_STAT) LBL_MAX, grp_lo(ig), grp_hi(ig), &
            DASH14, DASH14, max2_grp, DASH14, max3_grp
        call OutLine(trim(ln))

        sum2_all = sum2_all + sum2_grp;  max2_all = max(max2_all, max2_grp)
        sum3_all = sum3_all + sum3_grp;  max3_all = max(max3_all, max3_grp)

        if(ig < ngroup) call OutLine('  '//repeat('-',80))
    end do

    if(present(ret_avg_rd_list)) ret_avg_rd_list = sum2_all / ntot
    if(present(ret_max_rd_list)) ret_max_rd_list = max2_all
    if(present(ret_avg_rd_part)) ret_avg_rd_part = sum3_all / ntot
    if(present(ret_max_rd_part)) ret_max_rd_part = max3_all
end subroutine PrintTableGen

!  Total wall-clock timing summary across all cases.
subroutine PrintTotalSummary()
    real(dp) :: tot_1pt, tot_list, tot_part
    character(len=120) :: ln

    tot_1pt  = tF_1pt  + tC_1pt  + tE_1pt
    tot_list = tF_list + tC_list + tE_list
    tot_part = tF_part + tC_part + tE_part

    call OutLine('')
    call OutLine(repeat('=',100))
    call OutLine('  TOTAL TIMING SUMMARY (wall clock)')
    call OutLine('  Bins evaluated: '// &
        'A-fine ('//trim(i2s(NF))//') + '// &
        'A-coarse ('//trim(i2s(NC))//') + '// &
        'B ('//trim(i2s(NE))//') = '// &
        trim(i2s(NF+NC+NE))//' total')
    call OutLine(repeat('-',100))
    call OutLine('  Method               A-fine [s]   A-coarse [s]   B [s]        TOTAL [s]')
    call OutLine('  '//repeat('-',80))
    write(ln,'(2X,A,T23,4(F10.3,3X))') &
        'single-point loop', tF_1pt, tC_1pt, tE_1pt, tot_1pt
    call OutLine(trim(ln))
    write(ln,'(2X,A,T23,4(F10.3,3X))') &
        'list (no partition)', tF_list, tC_list, tE_list, tot_list
    call OutLine(trim(ln))
    write(ln,'(2X,A,T23,4(F10.3,3X))') &
        'list (partitioned)', tF_part, tC_part, tE_part, tot_part
    call OutLine(trim(ln))
    call OutLine('  '//repeat('-',80))
    write(ln,'(2X,A,T50,F7.2)') &
        'Overall speedup  list / 1pt:', tot_1pt / max(tot_list, 1.e-9_dp)
    call OutLine(trim(ln))
    write(ln,'(2X,A,T50,F7.2)') &
        'Overall speedup  part / 1pt:', tot_1pt / max(tot_part, 1.e-9_dp)
    call OutLine(trim(ln))
    call OutLine(repeat('=',100))
end subroutine PrintTotalSummary

!  Total precision summary across all cases.
!  rd = |X_method - X_1pt| / |X_1pt|
subroutine PrintTotalPrecisionSummary()
    character(len=120) :: ln
    real(dp) :: ovr_avg_list, ovr_avg_part

    ovr_avg_list = (prec_avg_list(1)*NF + prec_avg_list(2)*NC + prec_avg_list(3)*NE) &
                   / real(NF+NC+NE, dp)
    ovr_avg_part = (prec_avg_part(1)*NF + prec_avg_part(2)*NC + prec_avg_part(3)*NE) &
                   / real(NF+NC+NE, dp)

    call OutLine('')
    call OutLine(repeat('=',100))
    call OutLine('  TOTAL PRECISION SUMMARY')
    call OutLine('  rd = |X_method - X_1pt| / |X_1pt|  (X_1pt = single-point reference)')
    call OutLine(repeat('-',100))
    call OutLine('  Case                  list (no partition)              list (partitioned)')
    call OutLine('                         avg rd      max rd               avg rd      max rd')
    call OutLine('  '//repeat('-',80))
    write(ln,'(2X,A,T24,ES10.3,2X,ES10.3,T56,ES10.3,2X,ES10.3)') &
        'A-fine  ('//trim(i2s(NF))//' bins)', &
        prec_avg_list(1), prec_max_list(1), prec_avg_part(1), prec_max_part(1)
    call OutLine(trim(ln))
    write(ln,'(2X,A,T24,ES10.3,2X,ES10.3,T56,ES10.3,2X,ES10.3)') &
        'A-coarse ('//trim(i2s(NC))//' bins)', &
        prec_avg_list(2), prec_max_list(2), prec_avg_part(2), prec_max_part(2)
    call OutLine(trim(ln))
    write(ln,'(2X,A,T24,ES10.3,2X,ES10.3,T56,ES10.3,2X,ES10.3)') &
        'B       ('//trim(i2s(NE))//' bins)', &
        prec_avg_list(3), prec_max_list(3), prec_avg_part(3), prec_max_part(3)
    call OutLine(trim(ln))
    call OutLine('  '//repeat('-',80))
    write(ln,'(2X,A,T24,ES10.3,2X,ES10.3,T56,ES10.3,2X,ES10.3)') &
        'OVERALL  ('//trim(i2s(NF+NC+NE))//' bins)', &
        ovr_avg_list, maxval(prec_max_list), ovr_avg_part, maxval(prec_max_part)
    call OutLine(trim(ln))
    call OutLine(repeat('=',100))
end subroutine PrintTotalPrecisionSummary

function i2s(n) result(s)
    integer, intent(in) :: n
    character(len=12) :: s
    write(s,'(I0)') n
end function i2s

end program test_TMDX_DY
