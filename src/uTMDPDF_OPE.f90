!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.00
!
!	Evaluation of the small-b OPE for uTMDPDF
!	
!	if you use this module please, quote ????.????
!
!	ver 3.00: release (AV, 18.07.2023)
!
!				A.Vladimirov (18.07.2023)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!---- Concept---
! This module
! * computes the tw2-convolution C x f[x,b] at the optimal point
! * creates and saves the grids for it
! * only global variables are kept here
! * the most part of the code is universal, and shared by many such modules

module uTMDPDF_OPE
use aTMDe_Numerics
use IntegrationRoutines
use IO_functions
use QCDinput
use uTMDPDF_model_new
implicit none

!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v3.00"
character (len=11),parameter :: moduleName="uTMDPDF_OPE"
!Last appropriate version of constants-file
integer,parameter::inputver=0

!--- general
logical:: started=.false.
!! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
integer::outputLevel=2
!! variable that count number of WARNING mesagges. In order not to spam too much
integer::messageTrigger=6
integer::messageCounter=0 !!! actual counter

!!!------------------------- PARAMETERS DEFINED IN THE INI-file--------------

!!! Perturbative order
integer :: orderMain=2 !! LO=0, NLO=1,...

!!!! X-Grid parameters
!! over x: i=0...Nx, x_0=xMin
real(dp) :: xMin=0.00001_dp !!! min x 
integer :: Nx=400 !!! number of points in grid
real(dp) :: DeltaX !!! increment of grid
real(dp) :: parX=2._dp !!! parameter of the grid

!!!! B-Grid parameters
!! over b: i=0...Nb, x_nB=BMax
real(dp) :: BMAX=25._dp !!! maximum B
real(dp) :: BMIN=1d-6 !!! minimum B
integer :: NB=250 !!! number of points in grid
real(dp) :: DeltaB !!! increment of grid

!!!! Numerical parameters
real(dp) :: toleranceINT=1d-6  !!! tolerance for numerical integration
integer :: maxIteration=4000   !!! maximum iteration in the integrals (not used at the moment)

logical(dp) :: IsMuYdependent = .true.  !!! if mu is y independent, computation is much(!) faster

!!!! grid preparation
logical :: useGrid=.true.  !!!idicator that grid must be prepared
logical :: withGluon=.false.   !!!indicator the gluon is needed in the grid
logical :: runTest=.false.   !!!trigger to run the test

!!!------------------------- HARD-CODED PARAMETERS ----------------------
!!! Coefficient lists
integer,parameter::parametrizationLength=37
!! { Log[1-x], log[1-x]^2, log[1-x]^3, log[1-x]^4, log[1-x]^5  !exact
!!   1/x, log[x]/x, Log[x]^2/x  !exact
!! Log[x], log[x]^2, Log[x]^3, log[x]^4, Log[x]^5 !exact
!! T0,...,T23 (Chebyshev polynomials) }
!! The Lmu^2 part is exact the later parts are fitted, but exact if posible (e.g. Lmu and Nf parts for q->q)
!!!!! TO DO: UPDATE TO EXACT VALUES [use 1809.07084]!!

!!!------------------------- DYNAMICAL-GLOBAL PARAMETERS -------------------
real(dp) :: c4_global=1_dp  !!! scale variation parameter
logical :: gridReady!!!!indicator that grid is ready to use. If it is .true., the TMD calculated from the grid

!!!------------------------- SPECIAL VARIABLES FOR GRID (used by TMDGrid-XB)------------------
real(dp), dimension(:,:,:,:), allocatable :: gridMain !!!! THIS IS HUGE(!) matrix for the grid
real(dp), dimension(:,:,:,:), allocatable :: interpolationParameters !!!! for b>bGrid_Max we
integer, dimension(:), allocatable:: hadronsInGRID  !!! table that saves the number of hadons into plain list
integer::numberOfHadrons=1				!!!total number of hadrons to be stored

!!--------------------------------------Public interface-----------------------------------------
public::uTMDPDF_OPE_Initialize

!!!!!!----FOR TEST
public::MakeGrid,ExtractFromGrid,CxF_compute,TestGrid

contains

!! Coefficient function
INCLUDE 'Code/uTMDPDF/coeffFunc-new2.f90'
!INCLUDE 'Code/uTMDPDF/coeffFunc-new.f90'

!! X-grid routines
INCLUDE 'Code/Twist2/Twist2Xgrid.f90'
!! B-grid routines
INCLUDE 'Code/Twist2/Twist2Bgrid.f90'
!! Mellin convolution routine
INCLUDE 'Code/Twist2/Twist2Convolution-new.f90'
!! Grid operation
INCLUDE 'Code/Grids/TMDGrid-XB.f90'


function uTMDPDF_OPE_IsInitialized()
    logical::uTMDPDF_OPE_IsInitialized
    uTMDPDF_OPE_IsInitialized=started
end function uTMDPDF_OPE_IsInitialized

!! Initialization of the package
subroutine uTMDPDF_OPE_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path,line
    logical::initRequired
    character(len=8)::order_global
    logical::bSTAR_lambdaDependent
    integer::i,FILEver

    if(started) return

    if(.not.QCDinput_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing QCDinput (from ',moduleName,')'
        if(present(prefix)) then
            call QCDinput_Initialize(file,prefix)
        else
            call QCDinput_Initialize(file)
        end if
    end if

    if(present(prefix)) then
        path=trim(adjustl(prefix))//trim(adjustr(file))
    else
        path=trim(adjustr(file))
    end if

    !----------------- reading ini-file --------------------------------------
    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    !!! Search for output level
    call MoveTO(51,'*0   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) FILEver
    if(FILEver<inputver) then
        write(*,*) 'artemide.'//trim(moduleName)//': const-file version is too old.'
        write(*,*) '		     Update the const-file with artemide.setup'
        write(*,*) '  '
        stop
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel
    if(outputLevel>1) write(*,*) '--------------------------------------------- '
    if(outputLevel>1) write(*,*) 'artemide.',moduleName," ",version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger

    !!!!!!!!!!! OPEN MP Initialization
    !$ call MoveTO(51,'*C   ')
    !$ call MoveTO(51,'*p1  ')
    !$read(51,*) i
    !$ call OMP_set_num_threads(i)

    call MoveTO(51,'*4   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>1) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        return
    end if

     !----- ORDER
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) order_global

    SELECT CASE(trim(order_global))
        CASE ("NA")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NA',color(" (TMD=fNP)",c_yellow)
            orderMain=-50
        CASE ("LO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: LO'
            orderMain=0
        CASE ("NLO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO'
            orderMain=1
        CASE ("NNLO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NNLO'
            orderMain=2
        CASE ("N2LO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NNLO'
            orderMain=2
        CASE ("NNNLO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: N3LO'
            orderMain=3
        CASE ("N3LO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: N3LO'
            orderMain=3
        CASE DEFAULT
            if(outputLevel>0)write(*,*) &
                WarningString('Initialize: unknown order for coefficient function. Switch to NLO.',moduleName)
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO'
            orderMain=1
        END SELECT

    if(outputLevel>2 .and. orderMain>-1) write(*,'(A,I1)') ' |  Coef.func.    =as^',orderMain

    call MoveTO(51,'*p2  ')
    read(51,*) useGrid
    call MoveTO(51,'*p3  ')
    read(51,*) withGluon
    call MoveTO(51,'*p4  ')
    read(51,*) numberOfHadrons
    call MoveTO(51,'*p5  ')
    read(51,*) runTest

    !!!!! ---- parameters of numerical evaluation
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) toleranceINT
    call MoveTO(51,'*p2  ')
    read(51,*) maxIteration

    if(outputLevel>2) then
        write(*,'(A,ES10.3)') ' |  tolerance     =',toleranceINT
        write(*,'(A,ES10.3)') ' |  max iteration =',REAL(maxIteration)
        end if

    !-------------Parameters of grid
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) xMin
    call MoveTO(51,'*p2  ')
    read(51,*) parX
    call MoveTO(51,'*p3  ')
    read(51,*) BMIN
    call MoveTO(51,'*p4  ')
    read(51,*) BMAX
    call MoveTO(51,'*p5  ')
    read(51,*) NX
    call MoveTO(51,'*p6  ')
    read(51,*) NB

    if(outputLevel>2) then
        write(*,*) 'Grid options:'
        write(*,'(A,ES10.3)') ' |  xGrid_Min                 =',xMin
        write(*,'(A,ES10.3,","ES10.3,"]")') ' |  bMIN        =[',bMIN,bMAX
        write(*,'(A,I6,A,I6,A)') ' |  (GridSizeX,GridSizeB)     =(',NX,',',NB,')'
        write(*,'(A,I3)')   ' |  hadrons to grid           =',numberOfHadrons
    end if

    CLOSE (51, STATUS='KEEP')

    c4_global=1d0

    !!! call initializations for Grids
    call XGrid_Initialize()
    call BGrid_Initialize()
    call TMDGrid_XB_Initialize()
    
    !!! Model initialisation is called from the uTMDPDF-module
    
    !!!!!!!Checking the x-dependance of muOPE
    IsMuYdependent=testMU()

    if(IsMuYdependent) then
        if(outputLevel>2) write(*,*) trim(moduleName)//': mu OPE is dependent on x'
    else
        if(outputLevel>2) write(*,*) trim(moduleName)//': mu OPE is independent on x'
    end if

    gridReady=.false.

    call MakeGrid()

    gridReady=.true.

    if(runTest) call TestGrid()

    started=.true.
    messageCounter=0

    if(outputLevel>0) write(*,*) color('----- arTeMiDe.uTMDPDF_OPE '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '

end subroutine uTMDPDF_OPE_Initialize

!!!!array of x times PDF(x,Q) for hadron 'hadron'
!!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function xf(x,Q,hadron)
    real(dp) :: x,Q
    integer:: hadron
    real(dp), dimension(-5:5):: xf
    
    xf=xPDF(x,Q,hadron)
    
end function xf

end module uTMDPDF_OPE
