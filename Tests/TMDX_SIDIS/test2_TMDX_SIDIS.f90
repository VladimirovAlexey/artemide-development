!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   test2_TMDX_SIDIS.f90
!
!   Smoothness test: verifies that the bin-averaged SIDIS cross-section converges to
!   the differential cross-section as each bin dimension is halved twice and then
!   collapsed to a single point (zero bin size).
!
!   Three process variants (all using the h1->h2 unpolarized SIDIS sub-process 2001):
!   Case A  proc1=1  d^5sigma/dQ^2 dz dx dpPerp^2    Q=[6,8] GeV, z=[0.3,0.5], x=[0.1,0.2]
!   Case B  proc1=2  d^5sigma/dy   dz dx dpPerp^2    y=[0.3,0.5] passed in Q-slot
!   Case C  proc1=3  d^5sigma/dQ^2 dz dy dpPerp^2    y=[0.3,0.5] passed in x-slot
!   All cases: pT=[0.3,0.5] GeV,  sqrt(s)=30 GeV,  no cuts
!
!   For each variable the bin shrinks in 4 steps:
!     step 0 -- full bin
!     step 1 -- halved  (bin_size/2, centered at reference)
!     step 2 -- quartered
!     step 3 -- point (lo = hi = reference value)
!
!   rd = |xSec(step) - xSec(point)| / |xSec(point)|
!
!   NP parameters: same as test_TMDX_SIDIS (artemide_SetNPparameters, combined uTMDPDF+uTMDFF)
!
!   Compile and run from the artemide root directory:
!     make program TARGET=Tests/TMDX_SIDIS/test2_TMDX_SIDIS.f90
!
!   The INI file must be present at Tests/TMDX_SIDIS/test.atmde.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test2_TMDX_SIDIS
use aTMDe_control
use TMDX_SIDIS
use aTMDe_Numerics
implicit none

character(*), parameter :: inifile = 'test.atmde'
character(*), parameter :: prefix  = 'Tests/TMDX_SIDIS/'
character(*), parameter :: outfile = 'Tests/TMDX_SIDIS/test2_TMDX_SIDIS.out'

real(dp), parameter :: NP_param(28) = [ &
    1.5_dp,     0.0859102_dp, 0.030294_dp,   0.0_dp,       &
    0.48622_dp, 0.0411753_dp, 0.569024_dp,   0.146933_dp,  &
    5.26034_dp, 21.1222_dp,   7.71185_dp,    0.1565_dp,    &
    0.240061_dp,0.0691505_dp, 1.0_dp,        1.0_dp,       &
    0.696096_dp,0.626588_dp,  0.00331303_dp,-0.466377_dp,  &
    0.88367_dp, 0.882092_dp,  1.74168_dp,    1.15036_dp,   &
    0.610318_dp,-0.101387_dp, 0.0_dp,        0.0_dp        ]

! Number of convergence steps: 0=full, 1=half, 2=quarter, 3=point
integer, parameter :: NS = 4

!=======================================================================
! Process codes and energy
!=======================================================================
integer,  parameter :: PROC_A(4) = [1, 1, 1, 2001]  ! d/dQ^2 dz dx dpPerp^2
integer,  parameter :: PROC_B(4) = [2, 1, 1, 2001]  ! d/dy   dz dx dpPerp^2  (Q-slot = y)
integer,  parameter :: PROC_C(4) = [3, 1, 1, 2001]  ! d/dQ^2 dz dy dpPerp^2  (x-slot = y)

real(dp), parameter :: s_SIDIS = 30._dp**2  ! 900 GeV^2
! Note: s=18^2=324 GeV^2 (as in test1) is too small for Q=[6,8] — gives y>1 at x=0.1.

real(dp), parameter :: NO_CUTS(4) = [0._dp, 0._dp, 0._dp, 0._dp]

!=======================================================================
! Reference bins (used as fixed ranges when a different variable is being tested)
!=======================================================================
real(dp), parameter :: Q_FIX(2)  = [6._dp,   8._dp ]   ! center Q_ref = 7 GeV
real(dp), parameter :: z_FIX(2)  = [0.3_dp,  0.5_dp]   ! center z_ref = 0.4
real(dp), parameter :: x_FIX(2)  = [0.1_dp,  0.2_dp]   ! center x_ref = 0.15
real(dp), parameter :: pT_FIX(2) = [0.3_dp,  0.5_dp]   ! center pT_ref = 0.4 GeV
real(dp), parameter :: y_FIX(2)  = [0.3_dp,  0.5_dp]   ! center y_ref  = 0.4 (proc1=2,3)

!=======================================================================
! Step arrays: each variable shrinks symmetrically toward its reference value.
!=======================================================================
! Q steps  -- center 7 GeV, half-width 1 GeV
real(dp), parameter :: Q_STEPS(NS,2) = reshape( &
    [6._dp,  6.5_dp, 6.75_dp, 7._dp, &
     8._dp,  7.5_dp, 7.25_dp, 7._dp], [NS,2])

! z steps  -- center 0.4, half-width 0.1
real(dp), parameter :: z_STEPS(NS,2) = reshape( &
    [0.3_dp, 0.35_dp, 0.375_dp, 0.4_dp, &
     0.5_dp, 0.45_dp, 0.425_dp, 0.4_dp], [NS,2])

! x steps  -- center 0.15, half-width 0.05
real(dp), parameter :: x_STEPS(NS,2) = reshape( &
    [0.1_dp,  0.125_dp, 0.1375_dp, 0.15_dp, &
     0.2_dp,  0.175_dp, 0.1625_dp, 0.15_dp], [NS,2])

! pT steps -- center 0.4 GeV, half-width 0.1 GeV
real(dp), parameter :: pT_STEPS(NS,2) = reshape( &
    [0.3_dp, 0.35_dp, 0.375_dp, 0.4_dp, &
     0.5_dp, 0.45_dp, 0.425_dp, 0.4_dp], [NS,2])

! y steps  -- center 0.4, half-width 0.1 (used for proc1=2 Q-slot and proc1=3 x-slot)
real(dp), parameter :: y_STEPS(NS,2) = reshape( &
    [0.3_dp, 0.35_dp, 0.375_dp, 0.4_dp, &
     0.5_dp, 0.45_dp, 0.425_dp, 0.4_dp], [NS,2])

!=======================================================================
integer :: fu, i
real(dp) :: vals(NS)

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
call PrintHeader()

!======================================================================
! CASE A: proc1=1   d^5sigma / dQ^2 dz dx dpPerp^2
!======================================================================
call SectionHeader('CASE A  proc1=1  process=[1,1,1,2001]  d^5s/dQ^2 dz dx dpPerp^2')

call SubHeader('Q-bin convergence  (fixed: z=[0.3,0.5]  x=[0.1,0.2]  pT=[0.3,0.5] GeV)')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_A, s_SIDIS, pT_FIX, z_FIX, x_FIX, Q_STEPS(i,:), .false., NO_CUTS)
end do
call PrintTable('Q', Q_STEPS, vals)

call SubHeader('z-bin convergence  (fixed: Q=[6,8] GeV  x=[0.1,0.2]  pT=[0.3,0.5] GeV)')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_A, s_SIDIS, pT_FIX, z_STEPS(i,:), x_FIX, Q_FIX, .false., NO_CUTS)
end do
call PrintTable('z', z_STEPS, vals)

call SubHeader('x-bin convergence  (fixed: Q=[6,8] GeV  z=[0.3,0.5]  pT=[0.3,0.5] GeV)')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_A, s_SIDIS, pT_FIX, z_FIX, x_STEPS(i,:), Q_FIX, .false., NO_CUTS)
end do
call PrintTable('x', x_STEPS, vals)

call SubHeader('pT-bin convergence  (fixed: Q=[6,8] GeV  z=[0.3,0.5]  x=[0.1,0.2])')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_A, s_SIDIS, pT_STEPS(i,:), z_FIX, x_FIX, Q_FIX, .false., NO_CUTS)
end do
call PrintTable('pT', pT_STEPS, vals)

call SubHeader('Simultaneous convergence  (Q, z, x, pT shrink together; ref: Q=7 z=0.4 x=0.15 pT=0.4 GeV)')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_A, s_SIDIS, pT_STEPS(i,:), z_STEPS(i,:), x_STEPS(i,:), Q_STEPS(i,:), &
                    .false., NO_CUTS)
end do
call PrintTableMulti4('Q', Q_STEPS, 'z', z_STEPS, 'x', x_STEPS, 'pT', pT_STEPS, vals)

!======================================================================
! CASE B: proc1=2   d^5sigma / dy dz dx dpPerp^2   (y passed in Q-slot)
!======================================================================
call SectionHeader('CASE B  proc1=2  process=[2,1,1,2001]  d^5s/dy dz dx dpPerp^2  (y in Q-slot)')

call SubHeader('y-bin convergence  (fixed: z=[0.3,0.5]  x=[0.1,0.2]  pT=[0.3,0.5] GeV)')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_B, s_SIDIS, pT_FIX, z_FIX, x_FIX, y_STEPS(i,:), .false., NO_CUTS)
end do
call PrintTable('y', y_STEPS, vals)

call SubHeader('z-bin convergence  (fixed: y=[0.3,0.5]  x=[0.1,0.2]  pT=[0.3,0.5] GeV)')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_B, s_SIDIS, pT_FIX, z_STEPS(i,:), x_FIX, y_FIX, .false., NO_CUTS)
end do
call PrintTable('z', z_STEPS, vals)

call SubHeader('x-bin convergence  (fixed: y=[0.3,0.5]  z=[0.3,0.5]  pT=[0.3,0.5] GeV)')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_B, s_SIDIS, pT_FIX, z_FIX, x_STEPS(i,:), y_FIX, .false., NO_CUTS)
end do
call PrintTable('x', x_STEPS, vals)

call SubHeader('pT-bin convergence  (fixed: y=[0.3,0.5]  z=[0.3,0.5]  x=[0.1,0.2])')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_B, s_SIDIS, pT_STEPS(i,:), z_FIX, x_FIX, y_FIX, .false., NO_CUTS)
end do
call PrintTable('pT', pT_STEPS, vals)

call SubHeader('Simultaneous convergence  (y, z, x, pT shrink together; ref: y=0.4 z=0.4 x=0.15 pT=0.4 GeV)')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_B, s_SIDIS, pT_STEPS(i,:), z_STEPS(i,:), x_STEPS(i,:), y_STEPS(i,:), &
                    .false., NO_CUTS)
end do
call PrintTableMulti4('y', y_STEPS, 'z', z_STEPS, 'x', x_STEPS, 'pT', pT_STEPS, vals)

!======================================================================
! CASE C: proc1=3   d^5sigma / dQ^2 dz dy dpPerp^2   (y passed in x-slot)
!======================================================================
call SectionHeader('CASE C  proc1=3  process=[3,1,1,2001]  d^5s/dQ^2 dz dy dpPerp^2  (y in x-slot)')

call SubHeader('Q-bin convergence  (fixed: z=[0.3,0.5]  y=[0.3,0.5]  pT=[0.3,0.5] GeV)')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_C, s_SIDIS, pT_FIX, z_FIX, y_FIX, Q_STEPS(i,:), .false., NO_CUTS)
end do
call PrintTable('Q', Q_STEPS, vals)

call SubHeader('z-bin convergence  (fixed: Q=[6,8] GeV  y=[0.3,0.5]  pT=[0.3,0.5] GeV)')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_C, s_SIDIS, pT_FIX, z_STEPS(i,:), y_FIX, Q_FIX, .false., NO_CUTS)
end do
call PrintTable('z', z_STEPS, vals)

call SubHeader('y-bin convergence  (fixed: Q=[6,8] GeV  z=[0.3,0.5]  pT=[0.3,0.5] GeV)')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_C, s_SIDIS, pT_FIX, z_FIX, y_STEPS(i,:), Q_FIX, .false., NO_CUTS)
end do
call PrintTable('y', y_STEPS, vals)

call SubHeader('pT-bin convergence  (fixed: Q=[6,8] GeV  z=[0.3,0.5]  y=[0.3,0.5])')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_C, s_SIDIS, pT_STEPS(i,:), z_FIX, y_FIX, Q_FIX, .false., NO_CUTS)
end do
call PrintTable('pT', pT_STEPS, vals)

call SubHeader('Simultaneous convergence  (Q, z, y, pT shrink together; ref: Q=7 z=0.4 y=0.4 pT=0.4 GeV)')
do i = 1, NS
    call xSec_SIDIS(vals(i), PROC_C, s_SIDIS, pT_STEPS(i,:), z_STEPS(i,:), y_STEPS(i,:), Q_STEPS(i,:), &
                    .false., NO_CUTS)
end do
call PrintTableMulti4('Q', Q_STEPS, 'z', z_STEPS, 'y', y_STEPS, 'pT', pT_STEPS, vals)

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
    call OutLine('  TMDX_SIDIS test2: smooth bin-to-point convergence')
    call OutLine('  Process: h1->h2 unpolarized SIDIS [?,1,1,2001],  sqrt(s)=30 GeV,  no cuts')
    call OutLine('  Case A  proc1=1  Q=[6,8] GeV  z=[0.3,0.5]  x=[0.1,0.2]  pT=[0.3,0.5] GeV')
    call OutLine('  Case B  proc1=2  y=[0.3,0.5] (Q-slot)  z=[0.3,0.5]  x=[0.1,0.2]  pT=[0.3,0.5] GeV')
    call OutLine('  Case C  proc1=3  Q=[6,8] GeV  z=[0.3,0.5]  y=[0.3,0.5] (x-slot)  pT=[0.3,0.5] GeV')
    call OutLine('  Each variable halved twice then collapsed to a point (step 3).')
    call OutLine('  rd = |xSec(step) - xSec(point)| / |xSec(point)|')
    call OutLine(repeat('=',100))
end subroutine PrintHeader

!  Prints a 4-row single-variable convergence table.
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

!  Prints a 4-row simultaneous-convergence table for 4 variables.
!
!  Uses F7.4 (7 chars) for the 8 variable columns (lo/hi pairs), with
!  1X gap between lo and hi of the same variable and 2X gap between variables.
!
!  Column layout (1-indexed character positions):
!    1-4   : 4X indent
!    5     : step digit
!    6-10  : 5X gap  (or "(pt) " for point row)
!    11-17 : V1_lo (F7.4)   18: 1X   19-25: V1_hi (F7.4)
!    26-27 : 2X gap
!    28-34 : V2_lo (F7.4)   35: 1X   36-42: V2_hi (F7.4)
!    43-44 : 2X gap
!    45-51 : V3_lo (F7.4)   52: 1X   53-59: V3_hi (F7.4)
!    60-61 : 2X gap
!    62-68 : V4_lo (F7.4)   69: 1X   70-76: V4_hi (F7.4)
!    77-78 : 2X gap
!    79-92 : xSec (ES14.6)
!    93-94 : 2X gap
!    95-103: rd (ES9.2)
subroutine PrintTableMulti4(v1n, v1s, v2n, v2s, v3n, v3s, v4n, v4s, vals)
    character(len=*), intent(in) :: v1n, v2n, v3n, v4n
    real(dp),         intent(in) :: v1s(NS,2), v2s(NS,2), v3s(NS,2), v4s(NS,2), vals(NS)
    integer  :: i
    real(dp) :: rd
    character(len=200) :: ln
    character(len=7)   :: l1, h1, l2, h2, l3, h3, l4, h4, rdname7
    character(len=14)  :: xname

    character(*), parameter :: FMT_HEAD = &
      '(4X,"step  ",A7,1X,A7,2X,A7,1X,A7,2X,A7,1X,A7,2X,A7,1X,A7,2X,A14,2X,A9)'
    character(*), parameter :: FMT_ROW  = &
      '(4X,I1,5X,F7.4,1X,F7.4,2X,F7.4,1X,F7.4,2X,F7.4,1X,F7.4,2X,F7.4,1X,F7.4,2X,ES14.6,2X,ES9.2)'
    character(*), parameter :: FMT_PT   = &
      '(4X,I1,"(pt) ",F7.4,1X,F7.4,2X,F7.4,1X,F7.4,2X,F7.4,1X,F7.4,2X,F7.4,1X,F7.4,2X,ES14.6,2X,"(ref)    ")'

    l1 = adjustr(trim(v1n)//'_lo')
    h1 = adjustr(trim(v1n)//'_hi')
    l2 = adjustr(trim(v2n)//'_lo')
    h2 = adjustr(trim(v2n)//'_hi')
    l3 = adjustr(trim(v3n)//'_lo')
    h3 = adjustr(trim(v3n)//'_hi')
    l4 = adjustr(trim(v4n)//'_lo')
    h4 = adjustr(trim(v4n)//'_hi')
    xname   = adjustr('xSec [pb]')
    rdname7 = adjustr('rd')

    write(ln, FMT_HEAD) l1, h1, l2, h2, l3, h3, l4, h4, xname, rdname7
    call OutLine(trim(ln))
    call OutLine('  '//repeat('-',101))

    do i = 1, NS-1
        rd = abs(vals(i) - vals(NS)) / max(abs(vals(NS)), 1.e-300_dp)
        write(ln, FMT_ROW) i-1, &
            v1s(i,1), v1s(i,2), v2s(i,1), v2s(i,2), &
            v3s(i,1), v3s(i,2), v4s(i,1), v4s(i,2), &
            vals(i), rd
        call OutLine(trim(ln))
    end do
    write(ln, FMT_PT) NS-1, &
        v1s(NS,1), v1s(NS,2), v2s(NS,1), v2s(NS,2), &
        v3s(NS,1), v3s(NS,2), v4s(NS,1), v4s(NS,2), &
        vals(NS)
    call OutLine(trim(ln))
end subroutine PrintTableMulti4

end program test2_TMDX_SIDIS
