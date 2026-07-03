!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   test_TMDX_SIDIS.f90
!
!   Benchmark and comparison of three pPerp-spectrum evaluation methods:
!     (1) xSec_SIDIS      -- single-bin call, sequential loop
!     (2) xSec_SIDIS_List -- list call, no pPerp-partitioning
!     (3) xSec_SIDIS_List -- list call, with pPerp-partitioning (Chebyshev interpolation)
!
!   Process: h1->h2 unpolarized SIDIS [1,1,1,2001], sqrt(s)=18 GeV, no cuts
!     Q bins  : [2,3], [4,5], [8,9] GeV
!     x bins  : [0.01,0.05], [0.1,0.2], [0.3,0.4]
!     z bins  : [0.2,0.3], [0.4,0.5], [0.6,0.7]
!     pPerp   : 0 to 0.25*<Q> in steps of 0.1 GeV
!                 Q=[2,3]:  6 bins (0 to 0.6 GeV)
!                 Q=[4,5]: 11 bins (0 to 1.1 GeV)
!                 Q=[8,9]: 21 bins (0 to 2.1 GeV)
!
!   NP parameters: provided by user -- replace placeholder arrays below
!
!   Compile and run from the artemide root directory:
!     make program TARGET=Tests/TMDX_SIDIS/test_TMDX_SIDIS.f90
!
!   The INI file must be present at Tests/TMDX_SIDIS/test.atmde.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_TMDX_SIDIS
use aTMDe_control
use TMDX_SIDIS
use aTMDe_Numerics
implicit none

character(*), parameter :: inifile = 'test.atmde'
character(*), parameter :: prefix  = 'Tests/TMDX_SIDIS/'
character(*), parameter :: outfile = 'Tests/TMDX_SIDIS/test_TMDX_SIDIS.out'

!---- NP parameters: REPLACE with values provided separately ----
real(dp), parameter :: NP_param(31) = [ 1.5, 0.0859102, 0.030294, 0.0,&
                                        0.48622, 0.0411753, 0.569024, 0.146933,&
                                         5.26034, 21.1222, 7.71185, 0.1565,&
                                          0.240061, 0.0691505, 1.0, 1.0,&
                                           0.696096, 0.626588, 0.00331303, -0.466377,&
                                            0.88367, 0.882092, 1.74168, 1.15036,&
                                             0.610318, -0.101387, 0.0, 0.0]

!=======================================================================
! Process and kinematics
!=======================================================================
integer,  parameter :: PROC(4) = [1, 1, 1, 2001]
real(dp), parameter :: s_SIDIS  = 18._dp**2   ! 324 GeV^2

integer,  parameter :: NX = 3
integer,  parameter :: NZ = 3
integer,  parameter :: NQ = 3

real(dp), parameter :: XBINS(NX,2) = reshape( &
    [0.01_dp, 0.1_dp, 0.3_dp,   0.05_dp, 0.2_dp, 0.4_dp], [NX,2])
real(dp), parameter :: ZBINS(NZ,2) = reshape( &
    [0.2_dp, 0.4_dp, 0.6_dp,   0.3_dp, 0.5_dp, 0.7_dp], [NZ,2])
real(dp), parameter :: QBINS(NQ,2) = reshape( &
    [2._dp, 4._dp, 8._dp,   3._dp, 5._dp, 9._dp], [NQ,2])

real(dp), parameter :: DPT      = 0.1_dp

! pPerp range = 0.25*<Q>, N_PT = floor(0.25*<Q>/0.1):
!   <Q>[2,3]=2.5  -> 0.625 -> 6 bins
!   <Q>[4,5]=4.5  -> 1.125 -> 11 bins
!   <Q>[8,9]=8.5  -> 2.125 -> 21 bins
integer, parameter :: N_PT1 = 6
integer, parameter :: N_PT2 = 11
integer, parameter :: N_PT3 = 21

integer, parameter :: N1 = NX * NZ * N_PT1   ! 54 bins
integer, parameter :: N2 = NX * NZ * N_PT2   ! 99 bins
integer, parameter :: N3 = NX * NZ * N_PT3   ! 189 bins

!=======================================================================
! Arrays for Q=[2,3] (case 1), 54 bins
!=======================================================================
real(dp) :: X1_1pt(N1), X1_list(N1), X1_part(N1)
real(dp) :: pT_1(N1,2), Q_1(N1,2), x_1(N1,2), z_1(N1,2), sArr_1(N1), cutP_1(N1,4)
integer  :: proc_1(N1,4)
logical  :: cuts_1(N1)

!=======================================================================
! Arrays for Q=[4,5] (case 2), 99 bins
!=======================================================================
real(dp) :: X2_1pt(N2), X2_list(N2), X2_part(N2)
real(dp) :: pT_2(N2,2), Q_2(N2,2), x_2(N2,2), z_2(N2,2), sArr_2(N2), cutP_2(N2,4)
integer  :: proc_2(N2,4)
logical  :: cuts_2(N2)

!=======================================================================
! Arrays for Q=[8,9] (case 3), 189 bins
!=======================================================================
real(dp) :: X3_1pt(N3), X3_list(N3), X3_part(N3)
real(dp) :: pT_3(N3,2), Q_3(N3,2), x_3(N3,2), z_3(N3,2), sArr_3(N3), cutP_3(N3,4)
integer  :: proc_3(N3,4)
logical  :: cuts_3(N3)

!=======================================================================
! Accumulated precision stats from each PrintTableSIDIS call (3 cases)
!=======================================================================
real(dp) :: prec_avg_list(3), prec_max_list(3)
real(dp) :: prec_avg_part(3), prec_max_part(3)

real(dp) :: t1_1pt, t1_list, t1_part
real(dp) :: t2_1pt, t2_list, t2_part
real(dp) :: t3_1pt, t3_list, t3_part

integer :: fu, ix, iz, ipt, idx

!======================================================================
! Initialize
!======================================================================
call artemide_Initialize(inifile, prefix=prefix)
call artemide_SetNPparameters(NP_param)

if(.not. TMDX_SIDIS_IsInitialized()) then
    write(*,*) 'ERROR: TMDX_SIDIS is not initialized. Enable it in the ini file.'
    error stop
end if

open(newunit=fu, file=outfile, status='replace', action='write')

!======================================================================
! Build list arrays.
! Loop order: ix outer -> iz middle -> ipt inner
! This places consecutive pPerp bins for the same (x,z,Q) together,
! which is the pattern the partitioning algorithm exploits.
!======================================================================
do ix = 1, NX
    do iz = 1, NZ
        do ipt = 1, N_PT1
            idx = ((ix-1)*NZ + (iz-1))*N_PT1 + ipt
            pT_1(idx,1)   = (ipt-1)*DPT
            pT_1(idx,2)   =  ipt   *DPT
            Q_1(idx,1)    = QBINS(1,1);  Q_1(idx,2)  = QBINS(1,2)
            x_1(idx,1)    = XBINS(ix,1); x_1(idx,2)  = XBINS(ix,2)
            z_1(idx,1)    = ZBINS(iz,1); z_1(idx,2)  = ZBINS(iz,2)
            sArr_1(idx)   = s_SIDIS
            proc_1(idx,:) = PROC
            cuts_1(idx)   = .false.;  cutP_1(idx,:) = 0._dp
        end do
    end do
end do

do ix = 1, NX
    do iz = 1, NZ
        do ipt = 1, N_PT2
            idx = ((ix-1)*NZ + (iz-1))*N_PT2 + ipt
            pT_2(idx,1)   = (ipt-1)*DPT
            pT_2(idx,2)   =  ipt   *DPT
            Q_2(idx,1)    = QBINS(2,1);  Q_2(idx,2)  = QBINS(2,2)
            x_2(idx,1)    = XBINS(ix,1); x_2(idx,2)  = XBINS(ix,2)
            z_2(idx,1)    = ZBINS(iz,1); z_2(idx,2)  = ZBINS(iz,2)
            sArr_2(idx)   = s_SIDIS
            proc_2(idx,:) = PROC
            cuts_2(idx)   = .false.;  cutP_2(idx,:) = 0._dp
        end do
    end do
end do

do ix = 1, NX
    do iz = 1, NZ
        do ipt = 1, N_PT3
            idx = ((ix-1)*NZ + (iz-1))*N_PT3 + ipt
            pT_3(idx,1)   = (ipt-1)*DPT
            pT_3(idx,2)   =  ipt   *DPT
            Q_3(idx,1)    = QBINS(3,1);  Q_3(idx,2)  = QBINS(3,2)
            x_3(idx,1)    = XBINS(ix,1); x_3(idx,2)  = XBINS(ix,2)
            z_3(idx,1)    = ZBINS(iz,1); z_3(idx,2)  = ZBINS(iz,2)
            sArr_3(idx)   = s_SIDIS
            proc_3(idx,:) = PROC
            cuts_3(idx)   = .false.;  cutP_3(idx,:) = 0._dp
        end do
    end do
end do

call PrintHeader()

!======================================================================
! CASE 1 -- Q=[2,3] GeV: 6 pPerp bins per (x,z) group, 54 bins total
!======================================================================
call SectionHeader('CASE 1 -- Q=[2,3] GeV: N_pPerp=6 per (x,z) group, 54 bins total')

call system_clock_start(t1_1pt)
do idx = 1, N1
    call xSec_SIDIS(X1_1pt(idx), proc_1(idx,:), sArr_1(idx), &
                    pT_1(idx,:), z_1(idx,:), x_1(idx,:), Q_1(idx,:), &
                    cuts_1(idx), cutP_1(idx,:))
end do
call system_clock_stop(t1_1pt)

call system_clock_start(t1_list)
call xSec_SIDIS_List(X1_list, proc_1, sArr_1, pT_1, z_1, x_1, Q_1, cuts_1, cutP_1, &
                     doPartitioning=.false.)
call system_clock_stop(t1_list)

call system_clock_start(t1_part)
call xSec_SIDIS_List(X1_part, proc_1, sArr_1, pT_1, z_1, x_1, Q_1, cuts_1, cutP_1, &
                     doPartitioning=.true.)
call system_clock_stop(t1_part)

call PrintTimings(t1_1pt, t1_list, t1_part)
call PrintTableSIDIS(NX, NZ, N_PT1, N1, DPT, &
                     XBINS(:,1), XBINS(:,2), ZBINS(:,1), ZBINS(:,2), &
                     X1_1pt, X1_list, X1_part, &
                     prec_avg_list(1), prec_max_list(1), prec_avg_part(1), prec_max_part(1))

!======================================================================
! CASE 2 -- Q=[4,5] GeV: 11 pPerp bins per (x,z) group, 99 bins total
!======================================================================
call SectionHeader('CASE 2 -- Q=[4,5] GeV: N_pPerp=11 per (x,z) group, 99 bins total')

call system_clock_start(t2_1pt)
do idx = 1, N2
    call xSec_SIDIS(X2_1pt(idx), proc_2(idx,:), sArr_2(idx), &
                    pT_2(idx,:), z_2(idx,:), x_2(idx,:), Q_2(idx,:), &
                    cuts_2(idx), cutP_2(idx,:))
end do
call system_clock_stop(t2_1pt)

call system_clock_start(t2_list)
call xSec_SIDIS_List(X2_list, proc_2, sArr_2, pT_2, z_2, x_2, Q_2, cuts_2, cutP_2, &
                     doPartitioning=.false.)
call system_clock_stop(t2_list)

call system_clock_start(t2_part)
call xSec_SIDIS_List(X2_part, proc_2, sArr_2, pT_2, z_2, x_2, Q_2, cuts_2, cutP_2, &
                     doPartitioning=.true.)
call system_clock_stop(t2_part)

call PrintTimings(t2_1pt, t2_list, t2_part)
call PrintTableSIDIS(NX, NZ, N_PT2, N2, DPT, &
                     XBINS(:,1), XBINS(:,2), ZBINS(:,1), ZBINS(:,2), &
                     X2_1pt, X2_list, X2_part, &
                     prec_avg_list(2), prec_max_list(2), prec_avg_part(2), prec_max_part(2))

!======================================================================
! CASE 3 -- Q=[8,9] GeV: 21 pPerp bins per (x,z) group, 189 bins total
!======================================================================
call SectionHeader('CASE 3 -- Q=[8,9] GeV: N_pPerp=21 per (x,z) group, 189 bins total')

call system_clock_start(t3_1pt)
do idx = 1, N3
    call xSec_SIDIS(X3_1pt(idx), proc_3(idx,:), sArr_3(idx), &
                    pT_3(idx,:), z_3(idx,:), x_3(idx,:), Q_3(idx,:), &
                    cuts_3(idx), cutP_3(idx,:))
end do
call system_clock_stop(t3_1pt)

call system_clock_start(t3_list)
call xSec_SIDIS_List(X3_list, proc_3, sArr_3, pT_3, z_3, x_3, Q_3, cuts_3, cutP_3, &
                     doPartitioning=.false.)
call system_clock_stop(t3_list)

call system_clock_start(t3_part)
call xSec_SIDIS_List(X3_part, proc_3, sArr_3, pT_3, z_3, x_3, Q_3, cuts_3, cutP_3, &
                     doPartitioning=.true.)
call system_clock_stop(t3_part)

call PrintTimings(t3_1pt, t3_list, t3_part)
call PrintTableSIDIS(NX, NZ, N_PT3, N3, DPT, &
                     XBINS(:,1), XBINS(:,2), ZBINS(:,1), ZBINS(:,2), &
                     X3_1pt, X3_list, X3_part, &
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
    call OutLine(repeat('=',110))
    call OutLine('  '//trim(title))
    call OutLine(repeat('=',110))
end subroutine SectionHeader

subroutine PrintHeader()
    call OutLine(repeat('=',110))
    call OutLine('  TMDX_SIDIS test')
    call OutLine('  Process: h1->h2 unpolarized SIDIS [1,1,1,2001],  sqrt(s)=18 GeV,  no cuts')
    call OutLine('  Q bins : [2,3], [4,5], [8,9] GeV')
    call OutLine('  x bins : [0.01,0.05], [0.1,0.2], [0.3,0.4]')
    call OutLine('  z bins : [0.2,0.3], [0.4,0.5], [0.6,0.7]')
    call OutLine('  pPerp  : 0 to 0.25*<Q> in steps of 0.1 GeV  (N_pT = 6/11/21 per case)')
    call OutLine(repeat('=',110))
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

!  Table printer for SIDIS.
!
!  Column layout:
!    1-2   : 2X indent
!    3-17  : pPerp bin "[F5.2,F5.2]  " (15 chars)
!    18-34 : x bin "[F6.3,F6.3]  " (17 chars)
!    35-45 : z bin "[F4.2,F4.2]" (11 chars)
!    46-48 : 3X gap
!    49-62 : 1-point value (ES14.6)
!    63-65 : 3X gap
!    66-79 : list value (ES14.6)
!    80-82 : 3X gap
!    83-91 : rd(list) (ES9.2)
!    92-94 : 3X gap
!    95-108: partitioned value (ES14.6)
!    109-111: 3X gap
!    112-120: rd(part) (ES9.2)
!
!  rd = |X_method - X_1pt| / |X_1pt|  (X_1pt = single-point reference)
subroutine PrintTableSIDIS(nxg, nzg, npt, ntot, dpt, &
                           xb_lo, xb_hi, zb_lo, zb_hi, &
                           X1, X2, X3, &
                           ret_avg_rd_list, ret_max_rd_list, &
                           ret_avg_rd_part, ret_max_rd_part)
    integer,  intent(in) :: nxg, nzg, npt, ntot
    real(dp), intent(in) :: dpt
    real(dp), intent(in) :: xb_lo(nxg), xb_hi(nxg), zb_lo(nzg), zb_hi(nzg)
    real(dp), intent(in) :: X1(ntot), X2(ntot), X3(ntot)
    real(dp), intent(out), optional :: ret_avg_rd_list, ret_max_rd_list
    real(dp), intent(out), optional :: ret_avg_rd_part, ret_max_rd_part

    integer  :: ixg, izg, ipt, idx, ng
    real(dp) :: rd2, rd3
    real(dp) :: sum2_grp, max2_grp, sum3_grp, max3_grp
    real(dp) :: sum2_all, max2_all, sum3_all, max3_all
    character(len=200) :: ln

    character(len=15), parameter :: LBL_AVG = '[avg rd]       '
    character(len=15), parameter :: LBL_MAX = '[max rd]       '
    character(len=14), parameter :: DASH14  = '--------------'

    character(*), parameter :: FMT_DATA = &
      '(2X,"[",F5.2,",",F5.2,"]  ","[",F6.3,",",F6.3,"]  ","[",F4.2,",",F4.2,"]",'// &
      '3X,ES14.6,3X,ES14.6,3X,ES9.2,3X,ES14.6,3X,ES9.2)'
    character(*), parameter :: FMT_STAT = &
      '(2X,A15,"[",F6.3,",",F6.3,"]  ","[",F4.2,",",F4.2,"]",'// &
      '3X,A14,3X,A14,3X,ES9.2,3X,A14,3X,ES9.2)'

    sum2_all = 0._dp;  max2_all = 0._dp
    sum3_all = 0._dp;  max3_all = 0._dp

    call OutLine('')
    call OutLine( &
      '  rd = |X_method - X_1pt| / |X_1pt|  (X_1pt = single-point sequential reference)')
    call OutLine( &
      '  pPerp [GeV]     x bin              z bin       '// &
      '  1-point [pb]      list [pb]    rd(list)'// &
      '    partitioned [pb]   rd(part)')
    call OutLine('  '//repeat('-',130))

    do ixg = 1, nxg
        do izg = 1, nzg
            sum2_grp = 0._dp;  max2_grp = 0._dp
            sum3_grp = 0._dp;  max3_grp = 0._dp

            do ipt = 1, npt
                idx = ((ixg-1)*nzg + (izg-1))*npt + ipt
                rd2 = abs(X2(idx)-X1(idx)) / max(abs(X1(idx)), 1.e-300_dp)
                rd3 = abs(X3(idx)-X1(idx)) / max(abs(X1(idx)), 1.e-300_dp)
                write(ln, FMT_DATA) &
                    (ipt-1)*dpt, ipt*dpt, xb_lo(ixg), xb_hi(ixg), zb_lo(izg), zb_hi(izg), &
                    X1(idx), X2(idx), rd2, X3(idx), rd3
                call OutLine(trim(ln))
                sum2_grp = sum2_grp + rd2;  max2_grp = max(max2_grp, rd2)
                sum3_grp = sum3_grp + rd3;  max3_grp = max(max3_grp, rd3)
            end do

            write(ln, FMT_STAT) LBL_AVG, xb_lo(ixg), xb_hi(ixg), zb_lo(izg), zb_hi(izg), &
                DASH14, DASH14, sum2_grp/npt, DASH14, sum3_grp/npt
            call OutLine(trim(ln))
            write(ln, FMT_STAT) LBL_MAX, xb_lo(ixg), xb_hi(ixg), zb_lo(izg), zb_hi(izg), &
                DASH14, DASH14, max2_grp, DASH14, max3_grp
            call OutLine(trim(ln))

            sum2_all = sum2_all + sum2_grp;  max2_all = max(max2_all, max2_grp)
            sum3_all = sum3_all + sum3_grp;  max3_all = max(max3_all, max3_grp)

            ng = (ixg-1)*nzg + izg
            if(ng < nxg*nzg) call OutLine('  '//repeat('-',80))
        end do
    end do

    if(present(ret_avg_rd_list)) ret_avg_rd_list = sum2_all / ntot
    if(present(ret_max_rd_list)) ret_max_rd_list = max2_all
    if(present(ret_avg_rd_part)) ret_avg_rd_part = sum3_all / ntot
    if(present(ret_max_rd_part)) ret_max_rd_part = max3_all
end subroutine PrintTableSIDIS

subroutine PrintTotalSummary()
    real(dp) :: tot_1pt, tot_list, tot_part
    character(len=120) :: ln

    tot_1pt  = t1_1pt  + t2_1pt  + t3_1pt
    tot_list = t1_list + t2_list + t3_list
    tot_part = t1_part + t2_part + t3_part

    call OutLine('')
    call OutLine(repeat('=',110))
    call OutLine('  TOTAL TIMING SUMMARY (wall clock)')
    call OutLine('  Bins evaluated: '// &
        'Q=[2,3] ('//trim(i2s(N1))//') + '// &
        'Q=[4,5] ('//trim(i2s(N2))//') + '// &
        'Q=[8,9] ('//trim(i2s(N3))//') = '// &
        trim(i2s(N1+N2+N3))//' total')
    call OutLine(repeat('-',110))
    call OutLine('  Method               Q=[2,3] [s]    Q=[4,5] [s]    Q=[8,9] [s]    TOTAL [s]')
    call OutLine('  '//repeat('-',90))
    write(ln,'(2X,A,T23,4(F10.3,3X))') &
        'single-point loop', t1_1pt, t2_1pt, t3_1pt, tot_1pt
    call OutLine(trim(ln))
    write(ln,'(2X,A,T23,4(F10.3,3X))') &
        'list (no partition)', t1_list, t2_list, t3_list, tot_list
    call OutLine(trim(ln))
    write(ln,'(2X,A,T23,4(F10.3,3X))') &
        'list (partitioned)', t1_part, t2_part, t3_part, tot_part
    call OutLine(trim(ln))
    call OutLine('  '//repeat('-',90))
    write(ln,'(2X,A,T55,F7.2)') &
        'Overall speedup  list / 1pt:', tot_1pt / max(tot_list, 1.e-9_dp)
    call OutLine(trim(ln))
    write(ln,'(2X,A,T55,F7.2)') &
        'Overall speedup  part / 1pt:', tot_1pt / max(tot_part, 1.e-9_dp)
    call OutLine(trim(ln))
    call OutLine(repeat('=',110))
end subroutine PrintTotalSummary

subroutine PrintTotalPrecisionSummary()
    character(len=120) :: ln
    real(dp) :: ovr_avg_list, ovr_avg_part

    ovr_avg_list = (prec_avg_list(1)*N1 + prec_avg_list(2)*N2 + prec_avg_list(3)*N3) &
                   / real(N1+N2+N3, dp)
    ovr_avg_part = (prec_avg_part(1)*N1 + prec_avg_part(2)*N2 + prec_avg_part(3)*N3) &
                   / real(N1+N2+N3, dp)

    call OutLine('')
    call OutLine(repeat('=',110))
    call OutLine('  TOTAL PRECISION SUMMARY')
    call OutLine('  rd = |X_method - X_1pt| / |X_1pt|  (X_1pt = single-point reference)')
    call OutLine(repeat('-',110))
    call OutLine('  Case                  list (no partition)              list (partitioned)')
    call OutLine('                         avg rd      max rd               avg rd      max rd')
    call OutLine('  '//repeat('-',90))
    write(ln,'(2X,A,T24,ES10.3,2X,ES10.3,T56,ES10.3,2X,ES10.3)') &
        'Q=[2,3] ('//trim(i2s(N1))//' bins)', &
        prec_avg_list(1), prec_max_list(1), prec_avg_part(1), prec_max_part(1)
    call OutLine(trim(ln))
    write(ln,'(2X,A,T24,ES10.3,2X,ES10.3,T56,ES10.3,2X,ES10.3)') &
        'Q=[4,5] ('//trim(i2s(N2))//' bins)', &
        prec_avg_list(2), prec_max_list(2), prec_avg_part(2), prec_max_part(2)
    call OutLine(trim(ln))
    write(ln,'(2X,A,T24,ES10.3,2X,ES10.3,T56,ES10.3,2X,ES10.3)') &
        'Q=[8,9] ('//trim(i2s(N3))//' bins)', &
        prec_avg_list(3), prec_max_list(3), prec_avg_part(3), prec_max_part(3)
    call OutLine(trim(ln))
    call OutLine('  '//repeat('-',90))
    write(ln,'(2X,A,T24,ES10.3,2X,ES10.3,T56,ES10.3,2X,ES10.3)') &
        'OVERALL  ('//trim(i2s(N1+N2+N3))//' bins)', &
        ovr_avg_list, maxval(prec_max_list), ovr_avg_part, maxval(prec_max_part)
    call OutLine(trim(ln))
    call OutLine(repeat('=',110))
end subroutine PrintTotalPrecisionSummary

function i2s(n) result(s)
    integer, intent(in) :: n
    character(len=12) :: s
    write(s,'(I0)') n
end function i2s

end program test_TMDX_SIDIS
