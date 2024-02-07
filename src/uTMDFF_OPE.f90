!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.00
!
!    Evaluation of the small-b OPE for uTMDFF
!
!    if you use this module please, quote ????.????
!
!    ver 3.00: release (AV, 20.07.2023)
!
!                A.Vladimirov (20.07.2023)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!---- Concept---
! This module
! * computes the tw2-convolution C x f[x,b] at the optimal point
! * creates and saves the grids for it
! * only global variables are kept here
! * the most part of the code is universal, and shared by many such modules

module uTMDFF_OPE
use aTMDe_Numerics
use IntegrationRoutines
use IO_functions
use QCDinput
use uTMDFF_model
implicit none

!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v3.00"
character (len=11),parameter :: moduleName="uTMDFF_OPE"
!Last appropriate version of constants-file
integer,parameter::inputver=30

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
real(dp) :: toleranceGEN=1d-6  !!! tolerance for other purposes
integer :: maxIteration=4000   !!! maximum iteration in the integrals (not used at the moment)

logical(dp) :: IsMuYdependent = .true.  !!! if mu is y independent, computation is much(!) faster

!!!! grid preparation
logical :: useGrid=.true.  !!!idicator that grid must be prepared
logical :: withGluon=.false.   !!!indicator the gluon is needed in the grid
logical :: runTest=.false.   !!!trigger to run the test

!!!------------------------- HARD-CODED PARAMETERS ----------------------
!!! Coefficient lists
integer,parameter::parametrizationLength=36

!!!------------------------- DYNAMICAL-GLOBAL PARAMETERS -------------------
real(dp) :: c4_global=1_dp  !!! scale variation parameter
logical :: gridReady!!!!indicator that grid is ready to use. If it is .true., the TMD calculated from the grid

!!!------------------------- SPECIAL VARIABLES FOR GRID (used by TMDGrid-XB)------------------
real(dp), dimension(:,:,:,:), allocatable :: gridMain !!!! THIS IS HUGE(!) matrix for the grid
real(dp), dimension(:,:,:,:), allocatable :: interpolationParameters !!!! for b>bGrid_Max we interpolate
integer::numberOfHadrons=1                !!!total number of hadrons to be stored

!!--------------------------------------Public interface-----------------------------------------
public::uTMDFF_OPE_IsInitialized,uTMDFF_OPE_Initialize,uTMDFF_OPE_convolution
public::uTMDFF_OPE_resetGrid,uTMDFF_OPE_testGrid,uTMDFF_OPE_SetPDFreplica,uTMDFF_OPE_SetScaleVariation

!!!!!!----FOR TEST
!public::MakeGrid,ExtractFromGrid,CxF_compute,TestGrid

contains

!! Coefficient function
INCLUDE 'Code/uTMDFF/coeffFunc.f90'

!! X-grid routines
INCLUDE 'Code/Twist2/Twist2Xgrid.f90'
!! B-grid routines
INCLUDE 'Code/Twist2/Twist2Bgrid.f90'
!! Mellin convolution routine
INCLUDE 'Code/Twist2/Twist2Convolution.f90'
!! Grid operation
INCLUDE 'Code/Twist2/Twist2Grid-XB.f90'


function uTMDFF_OPE_IsInitialized()
    logical::uTMDFF_OPE_IsInitialized
    uTMDFF_OPE_IsInitialized=started
end function uTMDFF_OPE_IsInitialized

!! Initialization of the package
subroutine uTMDFF_OPE_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path
    logical::initRequired
    character(len=8)::order_global
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
        write(*,*) '             Update the const-file with artemide.setup'
        write(*,*) '  '
        CLOSE (51, STATUS='KEEP')
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

    call MoveTO(51,'*5   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>1) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        CLOSE (51, STATUS='KEEP')
        return
    end if

    !----- GENERALS
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) withGluon
    call MoveTO(51,'*p2  ')
    read(51,*) numberOfHadrons

    !----- ORDER
    call MoveTO(51,'*B   ')
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
    read(51,*) runTest

    !!!!! ---- parameters of numerical evaluation
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) toleranceINT
    call MoveTO(51,'*p2  ')
    read(51,*) toleranceGEN
    call MoveTO(51,'*p3  ')
    read(51,*) maxIteration

    if(outputLevel>2) then
        write(*,'(A,ES10.3)') ' |  tolerance     =',toleranceINT
        write(*,'(A,ES10.3)') ' |  max iteration =',REAL(maxIteration)
        end if

    !-------------Parameters of grid
    call MoveTO(51,'*E   ')
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
    
    !!! Model initialisation is called from the uTMDFF-module
    
    !!!!!!!Checking the x-dependance of muOPE
    IsMuYdependent=testMU()

    if(IsMuYdependent) then
        if(outputLevel>2) write(*,*) trim(moduleName)//': mu OPE is dependent on x'
    else
        if(outputLevel>2) write(*,*) trim(moduleName)//': mu OPE is independent on x'
    end if
    gridReady=.false.

    if(useGrid) then
        call MakeGrid()
        gridReady=.true.

        if(runTest) call TestGrid()
    end if

    started=.true.
    messageCounter=0

    if(outputLevel>0) write(*,*) color('----- arTeMiDe.uTMDFF_OPE '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '

end subroutine uTMDFF_OPE_Initialize

!!!!!!!--------------------------- DEFINING ROUTINES ------------------------------------------

!!!!array of x times PDF(x,Q) for hadron 'hadron'
!!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function xf(x,Q,hadron)
    real(dp) :: x,Q
    integer:: hadron
    real(dp), dimension(-5:5):: xf
    
    xf=xFF(x,Q,hadron)
    
end function xf


!!!!! in the case of TMDFF one computes the convolution as
!!! int_z^1 dy/y C[z/y] d[y]/y^2 = 1/z^3 int_z^1 dy [y^2C](y) D[z/y],
!!! where D(z)=z d(z), and [y^2C](y) is tabulated in the coeff.functions.
!!! the Convolution integral computes int_z^1 dy [y^2C](y) D[z/y]. Thus one needs to devide by z^3
function uTMDFF_OPE_convolution(x,b,h,addGluon)
    real(dp),dimension(-5:5)::uTMDFF_OPE_convolution
    real(dp),intent(in)::x,b
    integer,intent(in)::h
    logical,optional,intent(in)::addGluon

    logical::gluon

    !!! check gluonity
    if(present(addGluon)) then
        gluon=addGluon
    else
        gluon=withGluon
    end if

    !!!! test for boundaries is done in uTMDFF_lowScale5 (on the enty to this procedure)

    !!!! case NA
    if(orderMain==-50) then
        if(gluon) then
            uTMDFF_OPE_convolution=1._dp
        else
            uTMDFF_OPE_convolution=(/1._dp,1._dp,1._dp,1._dp,1._dp,0._dp,1._dp,1._dp,1._dp,1._dp,1._dp/)
        end if
        return
    end if

    !!! computation
    if(useGrid) then
        if(gridReady) then
            uTMDFF_OPE_convolution=ExtractFromGrid(x,b,h)/x**3
        else
            call Warning_Raise('Called OPE_convolution while grid is not ready.',messageCounter,messageTrigger,moduleName)
            call uTMDFF_OPE_resetGrid()
            uTMDFF_OPE_convolution=ExtractFromGrid(x,b,h)/x**3
        end if
    else
        uTMDFF_OPE_convolution=CxF_compute(x,b,h,gluon)/x**3
    end if

end function uTMDFF_OPE_convolution


!!!!!!!!!! ------------------------ SUPPORINTG ROUTINES --------------------------------------
!!! This subroutine force reconstruction of the grid (if griding is ON)
subroutine uTMDFF_OPE_resetGrid()
    gridReady=.false.
    if(useGrid) then
        if(outputLevel>1) write(*,*) 'arTeMiDe ',moduleName,':  Grid Reset. with c4=',c4_global
        call MakeGrid()

        gridReady=.true.
    end if
end subroutine uTMDFF_OPE_resetGrid

!!! This subroutine force reconstruction of the grid (if griding is ON)
subroutine uTMDFF_OPE_testGrid()
    if(useGrid) then
        if(.not.gridReady) then
            call uTMDFF_OPE_resetGrid()
        end if
        call TestGrid()
    end if
end subroutine uTMDFF_OPE_testGrid

!! call QCDinput to change the PDF replica number
!! unset the grid, since it should be recalculated fro different PDF replica.
subroutine uTMDFF_OPE_SetPDFreplica(rep,hadron)
    integer,intent(in):: rep,hadron
    logical::newPDF

    call QCDinput_SetPDFreplica(rep,hadron,newPDF)
    if(newPDF) then
        gridReady=.false.
        call uTMDFF_OPE_resetGrid()
    else
        if(outputLevel>1) write(*,"('arTeMiDe ',A,':  replica of FF (',I4,' is the same as the used one. Nothing is done!')") &
        moduleName, rep
    end if

end subroutine uTMDFF_OPE_SetPDFreplica

!!!! this routine set the variations of scales
!!!! it is used for the estimation of errors
subroutine uTMDFF_OPE_SetScaleVariation(c4_in)
    real(dp),intent(in)::c4_in
    if(c4_in<0.1d0 .or. c4_in>10.d0) then
        if(outputLevel>0) write(*,*) WarningString('variation in c4 is enourmous. c4 is set to 2',moduleName)
        c4_global=2d0
        call uTMDFF_OPE_resetGrid()
    else if(abs(c4_in-c4_global)<toleranceGEN) then
        if(outputLevel>1) write(*,*) color('uTMDFF: c4-variation is ignored. c4='//real8ToStr(c4_global),c_yellow)
    else
        c4_global=c4_in
        if(outputLevel>1) write(*,*) color('uTMDFF: set scale variations c4 as:'//real8ToStr(c4_global),c_yellow)
        call uTMDFF_OPE_resetGrid()
    end if
end subroutine uTMDFF_OPE_SetScaleVariation

end module uTMDFF_OPE
