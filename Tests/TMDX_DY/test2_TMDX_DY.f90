!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   test2_TMDX_DY.f90
!
!   Smoothness test: verifies that the bin-averaged cross-section converges to the
!   differential cross-section as the bin size is halved twice and then collapsed
!   to a single point (zero bin size).
!
!   Case A  proc1=1  pp->gamma*      sqrt(s)=39 GeV    process=[1,1,1,1]
!   Case B  proc1=2  pp->gamma*(xF)  sqrt(s)=39 GeV    process=[2,1,1,1]
!   Case C  proc1=1  pp->Z+gamma     sqrt(s)=13 TeV    process=[1,1,1,3]
!
!   For each variable the bin is:
!     step 0 -- full bin
!     step 1 -- halved  (bin_size / 2, centered at reference)
!     step 2 -- quartered (bin_size / 4)
!     step 3 -- point    (lo = hi = reference value)
!
!   rd = |xSec(step) - xSec(point)| / |xSec(point)|
!
!   NP parameters: uTMDPDF SV19 replica 0 (same as test_TMDX_DY)
!
!   Compile and run from the artemide root directory:
!     make program TARGET=Tests/TMDX_DY/test2_TMDX_DY.f90
!
!   The INI file must be present at Tests/TMDX_DY/test.atmde.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test2_TMDX_DY
use aTMDe_control
use TMDX_DY
use aTMDe_Numerics
implicit none

character(*), parameter :: inifile = 'test.atmde'
character(*), parameter :: prefix  = 'Tests/TMDX_DY/'
character(*), parameter :: outfile = 'Tests/TMDX_DY/test2_TMDX_DY.out'

real(dp), parameter :: NP(12) = [ &
    0.874245_dp, 0.913883_dp, 0.991563_dp, 6.05412_dp,  &
    0.353908_dp, 46.6064_dp,  0.115161_dp, 1.53235_dp,  &
    1.31966_dp,  0.434833_dp, 0.0_dp,      0.0_dp       ]

! Number of steps: 0=full, 1=half, 2=quarter, 3=point
integer, parameter :: NS = 4

!=======================================================================
! Case A/B: low-energy  pp->gamma*  sqrt(s)=39 GeV
!=======================================================================
real(dp), parameter :: s_E = 39._dp**2

integer,  parameter :: PROC1(4) = [1, 1, 1, 1]  ! proc1=1: rapidity (y)
integer,  parameter :: PROC2(4) = [2, 1, 1, 1]  ! proc1=2: xF

! Y steps (proc1=1), centered at y=0
real(dp), parameter :: Y_STEPS(NS,2) = reshape( &
    [-0.4_dp, -0.2_dp, -0.1_dp,  0._dp, &
      0.4_dp,  0.2_dp,  0.1_dp,  0._dp], [NS,2])

! xF steps (proc1=2), centered at xF=0.1
real(dp), parameter :: XF_STEPS(NS,2) = reshape( &
    [-0.1_dp,   0._dp, 0.05_dp, 0.1_dp, &
      0.3_dp,  0.2_dp, 0.15_dp, 0.1_dp], [NS,2])

! qT steps (cases A/B), centered at qT=0.5 GeV
real(dp), parameter :: QT_STEPS_E(NS,2) = reshape( &
    [0._dp,  0.25_dp, 0.375_dp, 0.5_dp, &
     1._dp,  0.75_dp, 0.625_dp, 0.5_dp], [NS,2])

! Q steps (cases A/B), centered at Q=11.5 GeV
real(dp), parameter :: Q_STEPS(NS,2) = reshape( &
    [11._dp, 11.25_dp, 11.375_dp, 11.5_dp, &
     12._dp, 11.75_dp, 11.625_dp, 11.5_dp], [NS,2])

! Fixed secondary ranges (cases A/B)
real(dp), parameter :: Q_FIX_E(2)  = [11._dp,  12._dp ]
real(dp), parameter :: QT_FIX_E(2) = [0.3_dp,  0.7_dp ]
real(dp), parameter :: Y_FIX_E(2)  = [-0.2_dp,  0.2_dp]
real(dp), parameter :: XF_FIX(2)   = [ 0._dp,   0.2_dp]

!=======================================================================
! Case C: Z-boson  pp->Z+gamma  sqrt(s)=13 TeV
!=======================================================================
real(dp), parameter :: s_Z = 13000._dp**2

integer,  parameter :: PROC_Z(4) = [1, 1, 1, 3]  ! proc1=1, Z+gamma

! Y steps (case C), centered at y=0
real(dp), parameter :: Y_STEPS_Z(NS,2) = reshape( &
    [-1._dp, -0.5_dp, -0.25_dp, 0._dp, &
      1._dp,  0.5_dp,  0.25_dp, 0._dp], [NS,2])

! qT steps (case C), centered at qT=5 GeV
real(dp), parameter :: QT_STEPS_Z(NS,2) = reshape( &
    [0._dp, 2.5_dp, 3.75_dp, 5._dp, &
    10._dp, 7.5_dp, 6.25_dp, 5._dp], [NS,2])

! Fixed secondary ranges (case C)
real(dp), parameter :: Q_FIX_Z(2)  = [70._dp, 110._dp]
real(dp), parameter :: QT_FIX_Z(2) = [ 3._dp,   7._dp]
real(dp), parameter :: Y_FIX_Z(2)  = [-0.5_dp,  0.5_dp]

!=======================================================================
! Simultaneous multi-variable step arrays
! All bins centered at their reference point; each step halves every dimension.
!=======================================================================
! Cases A/B: y/xF, qT, Q shrink simultaneously
!   y   reference = 0,    full half-width = 0.2
real(dp), parameter :: Y_MULTI_E(NS,2) = reshape( &
    [-0.2_dp, -0.1_dp, -0.05_dp, 0._dp, &
      0.2_dp,  0.1_dp,  0.05_dp, 0._dp], [NS,2])

!   xF  reference = 0.1,  full half-width = 0.1
real(dp), parameter :: XF_MULTI_E(NS,2) = reshape( &
    [0._dp,  0.05_dp, 0.075_dp, 0.1_dp, &
     0.2_dp, 0.15_dp, 0.125_dp, 0.1_dp], [NS,2])

!   qT  reference = 0.5,  full half-width = 0.2
real(dp), parameter :: QT_MULTI_E(NS,2) = reshape( &
    [0.3_dp, 0.4_dp, 0.45_dp, 0.5_dp, &
     0.7_dp, 0.6_dp, 0.55_dp, 0.5_dp], [NS,2])

!   Q   reference = 11.5, full half-width = 0.5
real(dp), parameter :: Q_MULTI_E(NS,2) = reshape( &
    [11._dp, 11.25_dp, 11.375_dp, 11.5_dp, &
     12._dp, 11.75_dp, 11.625_dp, 11.5_dp], [NS,2])

! Case C: y, qT shrink simultaneously; Q is kept fixed at [70, 110] GeV
!   y   reference = 0,    full half-width = 0.5
real(dp), parameter :: Y_MULTI_Z(NS,2) = reshape( &
    [-0.5_dp, -0.25_dp, -0.125_dp, 0._dp, &
      0.5_dp,  0.25_dp,  0.125_dp, 0._dp], [NS,2])

!   qT  reference = 5,    full half-width = 2
real(dp), parameter :: QT_MULTI_Z(NS,2) = reshape( &
    [3._dp, 4._dp, 4.5_dp, 5._dp, &
     7._dp, 6._dp, 5.5_dp, 5._dp], [NS,2])

!   Q   fixed at [70, 110] for all steps (note in table will clarify)
real(dp), parameter :: Q_MULTI_Z(NS,2) = reshape( &
    [ 70._dp,  70._dp,  70._dp,  70._dp, &
     110._dp, 110._dp, 110._dp, 110._dp], [NS,2])

!=======================================================================
integer :: fu, i
real(dp) :: vals(NS)

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
call PrintHeader()

!======================================================================
! CASE A: proc1=1  rapidity (y)
!======================================================================
call SectionHeader('CASE A  proc1=1  pp->gamma*  sqrt(s)=39 GeV  process=[1,1,1,1]')

call SubHeader('Y-bin convergence  (fixed: qT=[0.3,0.7] GeV, Q=[11,12] GeV)')
do i = 1, NS
    call xSec_DY(vals(i), PROC1, s_E, QT_FIX_E, Q_FIX_E, Y_STEPS(i,:), .false.)
end do
call PrintTable('y', Y_STEPS, vals)

call SubHeader('qT-bin convergence  (fixed: y=[-0.2,0.2], Q=[11,12] GeV)')
do i = 1, NS
    call xSec_DY(vals(i), PROC1, s_E, QT_STEPS_E(i,:), Q_FIX_E, Y_FIX_E, .false.)
end do
call PrintTable('qT', QT_STEPS_E, vals)

call SubHeader('Q-bin convergence  (fixed: y=[-0.2,0.2], qT=[0.3,0.7] GeV)')
do i = 1, NS
    call xSec_DY(vals(i), PROC1, s_E, QT_FIX_E, Q_STEPS(i,:), Y_FIX_E, .false.)
end do
call PrintTable('Q', Q_STEPS, vals)

call SubHeader('Simultaneous convergence  (y, qT, Q shrink together; ref: y=0, qT=0.5 GeV, Q=11.5 GeV)')
do i = 1, NS
    call xSec_DY(vals(i), PROC1, s_E, QT_MULTI_E(i,:), Q_MULTI_E(i,:), Y_MULTI_E(i,:), .false.)
end do
call PrintTableMulti('y', Y_MULTI_E, QT_MULTI_E, Q_MULTI_E, vals)

!======================================================================
! CASE B: proc1=2  xF
!======================================================================
call SectionHeader('CASE B  proc1=2  pp->gamma* (xF)  sqrt(s)=39 GeV  process=[2,1,1,1]')

call SubHeader('xF-bin convergence  (fixed: qT=[0.3,0.7] GeV, Q=[11,12] GeV)')
do i = 1, NS
    call xSec_DY(vals(i), PROC2, s_E, QT_FIX_E, Q_FIX_E, XF_STEPS(i,:), .false.)
end do
call PrintTable('xF', XF_STEPS, vals)

call SubHeader('qT-bin convergence  (fixed: xF=[0.0,0.2], Q=[11,12] GeV)')
do i = 1, NS
    call xSec_DY(vals(i), PROC2, s_E, QT_STEPS_E(i,:), Q_FIX_E, XF_FIX, .false.)
end do
call PrintTable('qT', QT_STEPS_E, vals)

call SubHeader('Q-bin convergence  (fixed: xF=[0.0,0.2], qT=[0.3,0.7] GeV)')
do i = 1, NS
    call xSec_DY(vals(i), PROC2, s_E, QT_FIX_E, Q_STEPS(i,:), XF_FIX, .false.)
end do
call PrintTable('Q', Q_STEPS, vals)

call SubHeader('Simultaneous convergence  (xF, qT, Q shrink together; ref: xF=0.1, qT=0.5 GeV, Q=11.5 GeV)')
do i = 1, NS
    call xSec_DY(vals(i), PROC2, s_E, QT_MULTI_E(i,:), Q_MULTI_E(i,:), XF_MULTI_E(i,:), .false.)
end do
call PrintTableMulti('xF', XF_MULTI_E, QT_MULTI_E, Q_MULTI_E, vals)

!======================================================================
! CASE C: Z-boson  proc1=1
!======================================================================
call SectionHeader('CASE C  proc1=1  pp->Z+gamma  sqrt(s)=13 TeV  process=[1,1,1,3]')

call SubHeader('Y-bin convergence  (fixed: qT=[3,7] GeV, Q=[70,110] GeV)')
do i = 1, NS
    call xSec_DY(vals(i), PROC_Z, s_Z, QT_FIX_Z, Q_FIX_Z, Y_STEPS_Z(i,:), .false.)
end do
call PrintTable('y', Y_STEPS_Z, vals)

call SubHeader('qT-bin convergence  (fixed: y=[-0.5,0.5], Q=[70,110] GeV)')
do i = 1, NS
    call xSec_DY(vals(i), PROC_Z, s_Z, QT_STEPS_Z(i,:), Q_FIX_Z, Y_FIX_Z, .false.)
end do
call PrintTable('qT', QT_STEPS_Z, vals)

call SubHeader('Simultaneous convergence  (y, qT shrink together; Q fixed at [70,110] GeV; ref: y=0, qT=5 GeV)')
do i = 1, NS
    call xSec_DY(vals(i), PROC_Z, s_Z, QT_MULTI_Z(i,:), Q_FIX_Z, Y_MULTI_Z(i,:), .false.)
end do
call PrintTableMulti('y', Y_MULTI_Z, QT_MULTI_Z, Q_MULTI_Z, vals)

close(fu)

!=======================================================================
contains
!=======================================================================

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

subroutine SubHeader(title)
    character(len=*), intent(in) :: title
    call OutLine('')
    call OutLine('  -- '//trim(title)//' --')
end subroutine SubHeader

subroutine PrintHeader()
    call OutLine(repeat('=',100))
    call OutLine('  TMDX_DY test2: smooth bin-to-point convergence')
    call OutLine('  NP (uTMDPDF SV19 rep.0):  0.874245  0.913883  0.991563  6.05412')
    call OutLine('                              0.353908  46.6064   0.115161  1.53235')
    call OutLine('                              1.31966   0.434833  0         0')
    call OutLine('  Case A  proc1=1  pp->gamma*    sqrt(s)=39 GeV   y ref=0,    qT ref=0.5 GeV, Q ref=11.5 GeV')
    call OutLine('  Case B  proc1=2  pp->gamma*    sqrt(s)=39 GeV   xF ref=0.1, qT ref=0.5 GeV, Q ref=11.5 GeV')
    call OutLine('  Case C  proc1=1  pp->Z+gamma   sqrt(s)=13 TeV   y ref=0,    qT ref=5 GeV,   Q in [70,110] GeV')
    call OutLine('  Each variable is halved twice then collapsed to a point (step 3).')
    call OutLine('  rd = |xSec(step) - xSec(point)| / |xSec(point)|')
    call OutLine(repeat('=',100))
end subroutine PrintHeader

!  Prints a 4-row convergence table with aligned columns.
!
!  Column layout (1-indexed character positions):
!    1-4   : 4X indent
!    5     : step digit
!    6-10  : 5X gap  (or "(pt) " for point row)
!    11-19 : lo value / label (F9.4 / A9 right-justified)
!    20-22 : 3X gap
!    23-31 : hi value / label (F9.4 / A9)
!    32-34 : 3X gap
!    35-43 : bin_size value / label (F9.4 / A9)
!    44-46 : 3X gap
!    47-60 : xSec value / label (ES14.6 / A14)
!    61-63 : 3X gap
!    64-72 : rd value / label (ES9.2 / A9)
subroutine PrintTable(varname, steps, vals)
    character(len=*), intent(in) :: varname
    real(dp),         intent(in) :: steps(NS,2), vals(NS)
    integer  :: i
    real(dp) :: rd, bsize
    character(len=200) :: ln
    character(len=9)   :: lname, hname, bname, rdname
    character(len=14)  :: xname

    character(*), parameter :: FMT_HEAD = &
      '(4X,"step  ",A9,3X,A9,3X,A9,3X,A14,3X,A9)'
    character(*), parameter :: FMT_ROW  = &
      '(4X,I1,5X,F9.4,3X,F9.4,3X,F9.4,3X,ES14.6,3X,ES9.2)'
    character(*), parameter :: FMT_PT   = &
      '(4X,I1,"(pt) ",F9.4,3X,F9.4,3X,F9.4,3X,ES14.6,3X,"(ref)    ")'

    ! Right-justify column labels within their field widths to match
    ! the right-justified numeric output of F9.4 / ES14.6 / ES9.2
    lname  = adjustr(trim(varname)//'_lo')
    hname  = adjustr(trim(varname)//'_hi')
    bname  = adjustr('bin_size')
    xname  = adjustr('xSec [pb]')
    rdname = adjustr('rd')

    write(ln, FMT_HEAD) lname, hname, bname, xname, rdname
    call OutLine(trim(ln))
    call OutLine('  '//repeat('-',88))

    do i = 1, NS-1
        bsize = steps(i,2) - steps(i,1)
        rd    = abs(vals(i) - vals(NS)) / max(abs(vals(NS)), 1.e-300_dp)
        write(ln, FMT_ROW) i-1, steps(i,1), steps(i,2), bsize, vals(i), rd
        call OutLine(trim(ln))
    end do
    bsize = steps(NS,2) - steps(NS,1)
    write(ln, FMT_PT) NS-1, steps(NS,1), steps(NS,2), bsize, vals(NS)
    call OutLine(trim(ln))
end subroutine PrintTable

!  Prints a 4-row convergence table with all three variables shown simultaneously.
!
!  Column layout (1-indexed character positions):
!    1-4   : 4X indent
!    5     : step digit
!    6-10  : 5X gap  (or "(pt) " for point row)
!    11-19 : v_lo    (F9.4 / A9)
!    20-21 : 2X gap
!    22-30 : v_hi    (F9.4 / A9)
!    31-32 : 2X gap
!    33-41 : qT_lo   (F9.4 / A9)
!    42-43 : 2X gap
!    44-52 : qT_hi   (F9.4 / A9)
!    53-54 : 2X gap
!    55-63 : Q_lo    (F9.4 / A9)
!    64-65 : 2X gap
!    66-74 : Q_hi    (F9.4 / A9)
!    75-76 : 2X gap
!    77-90 : xSec    (ES14.6 / A14)
!    91-92 : 2X gap
!    93-101: rd      (ES9.2 / A9)
subroutine PrintTableMulti(vname, v_steps, qt_steps, q_steps, vals)
    character(len=*), intent(in) :: vname
    real(dp),         intent(in) :: v_steps(NS,2), qt_steps(NS,2), q_steps(NS,2), vals(NS)
    integer  :: i
    real(dp) :: rd
    character(len=200) :: ln
    character(len=9)   :: lv, hv, lqt, hqt, lq, hq, rdname
    character(len=14)  :: xname

    character(*), parameter :: FMT_HEAD = &
      '(4X,"step  ",A9,2X,A9,2X,A9,2X,A9,2X,A9,2X,A9,2X,A14,2X,A9)'
    character(*), parameter :: FMT_ROW  = &
      '(4X,I1,5X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,ES14.6,2X,ES9.2)'
    character(*), parameter :: FMT_PT   = &
      '(4X,I1,"(pt) ",F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,F9.4,2X,ES14.6,2X,"(ref)    ")'

    lv     = adjustr(trim(vname)//'_lo')
    hv     = adjustr(trim(vname)//'_hi')
    lqt    = adjustr('qT_lo')
    hqt    = adjustr('qT_hi')
    lq     = adjustr('Q_lo')
    hq     = adjustr('Q_hi')
    xname  = adjustr('xSec [pb]')
    rdname = adjustr('rd')

    write(ln, FMT_HEAD) lv, hv, lqt, hqt, lq, hq, xname, rdname
    call OutLine(trim(ln))
    call OutLine('  '//repeat('-',99))

    do i = 1, NS-1
        rd = abs(vals(i) - vals(NS)) / max(abs(vals(NS)), 1.e-300_dp)
        write(ln, FMT_ROW) i-1, &
            v_steps(i,1),  v_steps(i,2), &
            qt_steps(i,1), qt_steps(i,2), &
            q_steps(i,1),  q_steps(i,2), &
            vals(i), rd
        call OutLine(trim(ln))
    end do
    write(ln, FMT_PT) NS-1, &
        v_steps(NS,1),  v_steps(NS,2), &
        qt_steps(NS,1), qt_steps(NS,2), &
        q_steps(NS,1),  q_steps(NS,2), &
        vals(NS)
    call OutLine(trim(ln))
end subroutine PrintTableMulti

end program test2_TMDX_DY
