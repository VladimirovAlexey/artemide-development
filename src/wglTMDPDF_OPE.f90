!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.04
!
!    Evaluation of the small-b OPE for wglTMDPDF
!
!    if you use this module please, quote 2209.00962
!
!    ver 3.00: release (AV, 07.01.2024)
!    ver 3.04: structural update to a global standard (AV, 25.06.2026)
!
!                A.Vladimirov (07.01.2024)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!---- Concept---
! This module
! * computes the tw2-convolution C x f[x,b] at the optimal point
! * creates and saves the grids for it
! * only global variables are kept here
! * the most part of the code is universal, and shared by many such modules

module wglTMDPDF_OPE
use aTMDe_Numerics
use aTMDe_Integration
use aTMDe_IO
use aTMDe_optGrid
use QCDinput
use TMD_AD, only : Dpert_atL
use wglTMDPDF_model
implicit none

!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v3.04"
character (len=13),parameter :: moduleName="wglTMDPDF_OPE"
!Last appropriate version of constants-file
integer,parameter::inputver=40

!--- general
logical:: started=.false.
!! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
!! 3=FULL output
integer::outputLevel=2
type(Warning_OBJ)::Warning_Handler

!!!------------------------- PARAMETERS DEFINED IN THE INI-file--------------

!!! Perturbative order
integer :: orderMainTW2=2 !! LO=0, NLO=1,...
integer :: orderMainTW3=-50 !! LO=0, NLO=1,...
!!! Order of large-X resummation
logical :: resumLargeX=.false. !!!! could not be resummation of largeX
integer :: orderLX=0 !! LO=0 [no-resummation], NLO=1,...

!!! Phase space limitations parameters
real(dp) :: xMin=0.001_dp !!! min x
real(dp) :: BMAX=25._dp !!! maximum B
real(dp) :: BMIN=1d-6 !!! minimum B

!!!! Numerical parameters
real(dp) :: toleranceINT=1d-6  !!! tolerance for numerical integration
real(dp) :: toleranceGEN=1d-6  !!! tolerance for other purposes
integer :: maxIteration=4000   !!! maximum iteration in the integrals (not used at the moment)

!!!total number of hadrons to be used
integer::numberOfHadrons=1

logical :: IsMuYdependent = .true.

!!!! grid preparation
logical :: useGridTW2=.true.  !!!indicator that grid for Tw2 part must be prepared
logical :: useGridTW3=.false.  !!!indicator that grid for Tw3 must be prepared
logical :: withGluon=.false.  !!!indicator that the gluon is needed in the grid (THIS IS COMMON PARAMETER FOR BOTH GRIDS)

type(optGrid)::mainGridTw2
type(optGrid)::mainGridTw3

!!!------------------------- HARD-CODED PARAMETERS ----------------------
!!! Coefficient lists
integer,parameter::parametrizationLength=4 !!![exact]

!!!------------------------- DYNAMICAL-GLOBAL PARAMETERS -------------------
real(dp) :: c4_tw2_global=1._dp  !!! scale variation parameter
real(dp) :: c4_tw3_global=1._dp  !!! scale variation parameter

!!--------------------------------------Public interface-----------------------------------------
public::wglTMDPDF_OPE_IsInitialized,wglTMDPDF_OPE_Initialize
public::wglTMDPDF_OPE_tw2,wglTMDPDF_OPE_tw2_resetGrid,wglTMDPDF_OPE_tw2_SetPDFreplica,wglTMDPDF_OPE_tw2_SetScaleVariation
public::wglTMDPDF_OPE_tw3,wglTMDPDF_OPE_tw3_resetGrid,wglTMDPDF_OPE_tw3_SetPDFreplica,wglTMDPDF_OPE_tw3_SetScaleVariation

!!!!!!----FOR TEST
!public::MakeGrid,ExtractFromGrid,CxF_compute,TestGrid

contains

!! Coefficient function for tw2 part plain
INCLUDE 'Code/wglTMDPDF/coeffFunc.f90'

!! Coefficient function for tw2 part large-X-resummed
INCLUDE 'Code/wglTMDPDF/coeffFunc_largeX.f90'

!! Mellin convolution routine
INCLUDE 'Code/Twist2/Twist2_WW.f90'


function wglTMDPDF_OPE_IsInitialized()
    logical::wglTMDPDF_OPE_IsInitialized
    wglTMDPDF_OPE_IsInitialized=started
end function wglTMDPDF_OPE_IsInitialized

!! Initialization of the package
subroutine wglTMDPDF_OPE_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=:),allocatable::path
    logical::initRequired
    character(len=8)::order_global
    integer::i,FILEver,messageTrigger
    real(dp),allocatable::subGridsX(:),subGridsB(:)

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
        path=trim(adjustl(prefix))//trim(adjustl(file))
    else
        path=trim(adjustl(file))
    end if

    !----------------- reading ini-file --------------------------------------
    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")

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

    call MoveTO(51,'*16  ')
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
            orderMainTW2=-50
        CASE ("LO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: LO'
            orderMainTW2=0
        CASE ("NLO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO'
            orderMainTW2=1
        CASE DEFAULT
            if(outputLevel>0)write(*,*) &
                WarningString('Initialize: unknown order for coefficient function. Switch to NLO.',moduleName)
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO'
            orderMainTW2=1
        END SELECT

    if(outputLevel>2 .and. orderMainTW2>-1) write(*,'(A,I1)') ' |  Coef.func.    =as^',orderMainTW2

    call MoveTO(51,'*p2  ')
    read(51,*) useGridTW2

    call MoveTO(51,'*p4  ')
    read(51,*) resumLargeX
    if(resumLargeX) then
        call MoveTO(51,'*p5  ')
        read(51,*) order_global

        SELECT CASE(trim(order_global))
            CASE ("LO")
                if(outputLevel>1) write(*,*) trim(moduleName)//' Large-X order set: LO'
                orderLX =0
            CASE ("NLO")
                if(outputLevel>1) write(*,*) trim(moduleName)//' Large-X order set: NLO'
                orderLX=1
            CASE ("NNLO")
                if(outputLevel>1) write(*,*) trim(moduleName)//' Large-X order set: NNLO'
                orderLX=2
            CASE ("N2LO")
                if(outputLevel>1) write(*,*) trim(moduleName)//' Large-X order set: NNLO'
                orderLX=2
            CASE ("NNNLO")
                if(outputLevel>1) write(*,*) trim(moduleName)//' Large-X order set: N3LO'
                orderLX=3
            CASE ("N3LO")
                if(outputLevel>1) write(*,*) trim(moduleName)//' Large-X order set: N3LO'
                orderLX=3
            CASE DEFAULT
                if(outputLevel>0) then
                    write(*,*) &
                    WarningString('Initialize: unknown order for large-X resummation of coefficient function.',moduleName)
                    write(*,*) WarningString('set to same order as the common part.',moduleName)
                end if
                orderLX=orderMainTW2
            END SELECT

    if(outputLevel>2 .and. orderLX>-1) write(*,'(A,I1)') ' |  Large-X       =as^',orderLX
    end if

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
    read(51,*) i
    allocate(subGridsX(0:i))
    call MoveTO(51,'*p2  ')
    read(51,*) subGridsX
    call MoveTO(51,'*p4  ')
    read(51,*) i
    allocate(subGridsB(0:i))
    call MoveTO(51,'*p5  ')
    read(51,*) subGridsB

    !-------------------tw3-part is placeholder for a moment
    call MoveTO(51,'*F   ')
    call MoveTO(51,'*p1  ')
    read(51,*) order_global

    SELECT CASE(trim(order_global))
        CASE ("NA")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order[tw3] set: NA',color(" (TMD=fNP)",c_yellow)
            orderMainTW3=-50
        CASE DEFAULT
            if(outputLevel>0)write(*,*) &
                WarningString('Initialize: the perturbative tw3 convolution is not implemented yet. Switch to NA.',moduleName)
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NA'
            orderMainTW3=-50
        END SELECT

    if(outputLevel>2 .and. orderMainTW3>-1) write(*,'(A,I1)') ' |  Coef.func.    =as^',orderMainTW3

    call MoveTO(51,'*p2  ')
    read(51,*) useGridTW3

    CLOSE (51, STATUS='KEEP')

    Warning_Handler=Warning_OBJ(moduleName=moduleName,messageCounter=0,messageTrigger=messageTrigger)

    c4_tw2_global=1d0
    c4_tw3_global=1d0


    xMin=subGridsX(0)
    bMin=subGridsB(0)
    bMax=subGridsB(size(subGridsB)-1)

    if(abs(subGridsX(size(subGridsX)-1)-1)>toleranceGEN) then
        error stop ErrorString("The last subgrid in X must complete by x=1. Initialization terminated",moduleName)
    end if

    mainGridTw2=optGrid(path,'*16  ','*E   ',numberOfHadrons,withGluon,moduleName,outputLevel)
    mainGridTw3=optGrid(path,'*16  ','*G   ',numberOfHadrons,withGluon,moduleName,outputLevel)
    
    !!! Model initialisation is called from the wglTMDPDF-module
    
    !!!!!!!Checking the x-dependence of muOPE
    IsMuYdependent=testMU()

    if(IsMuYdependent) then
        if(outputLevel>2) write(*,*) trim(moduleName)//': mu OPE is dependent on x'
    else
        if(outputLevel>2) write(*,*) trim(moduleName)//': mu OPE is independent on x'
    end if

    if(useGridTW2) then
        call mainGridTw2%MakeGrid(tw2_convolution)
    end if
    if(useGridTW3) then
        call mainGridTw3%MakeGrid(tw3_convolution)
    end if

    started=.true.

    if(outputLevel>0) write(*,*) color('----- arTeMiDe.wglTMDPDF_OPE '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '

end subroutine wglTMDPDF_OPE_Initialize

!!!!!! this is just interface of function CxF_compute and CxF_largeX_compute to optGrid
function tw2_convolution(x,bT,hadron)
    real(dp),dimension(-5:5)::tw2_convolution
    integer, intent(in)::hadron
    real(dp),intent(in)::x,bT

    !!!! case NA
    if(orderMainTW2==-50) then
        if(withGluon) then
            tw2_convolution=1._dp
        else
            tw2_convolution=(/1._dp,1._dp,1._dp,1._dp,1._dp,0._dp,1._dp,1._dp,1._dp,1._dp,1._dp/)
        end if
        return
    end if

    if(resumLargeX) then
        tw2_convolution=CxF_largeX_compute(x,bT,hadron,withGluon)
    else
        tw2_convolution=CxF_compute(x,bT,hadron,withGluon)
    end if
end function tw2_convolution

!!! This procedure computes the convolution of the tw3 PDF with the coefficient function
!!! It has interface suitable to store in the optGrid
!!! For now it is empty
function tw3_convolution(x,b,h)
    real(dp),dimension(-5:5)::tw3_convolution
    real(dp),intent(in)::x,b
    integer,intent(in)::h

    !!!! case NA
    if(orderMainTW3==-50) then
        if(withGluon) then
            tw3_convolution=1._dp
        else
            tw3_convolution=(/1._dp,1._dp,1._dp,1._dp,1._dp,0._dp,1._dp,1._dp,1._dp,1._dp,1._dp/)
        end if
        return
    end if

    !!!!! perturbative convolution is not implemented yet
    error stop ErrorString("Twist-3 convolution is not implemented yet!",moduleName)

end function tw3_convolution

!!!!!!!--------------------------- DEFINING ROUTINES ------------------------------------------

!!!!array of x times PDF(x,Q) for hadron 'hadron'
!!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function xf(x,Q,hadron)
    real(dp),intent(in) :: x,Q
    integer,intent(in):: hadron
    real(dp), dimension(-5:5):: xf
    
    xf=x_hPDF(x,Q,hadron)
    
end function xf


function wglTMDPDF_OPE_tw2(x,b,h)
    real(dp),dimension(-5:5)::wglTMDPDF_OPE_tw2
    real(dp),intent(in)::x,b
    integer,intent(in)::h

    !!!! The factor 1/x is replace by 1 because the WW part of wgl is
    !!!! -x^2\int dy/y^2 c[x/y]h[y] = x \int dy/y (-x/y c[x/y] h[y])
    !!!! so x is canceled, and the factor x/y is included in the coefficient function

    !!! computation
    if(useGridTW2) then
        wglTMDPDF_OPE_tw2=mainGridTw2%Extract(x,b,h)
    else
        wglTMDPDF_OPE_tw2=tw2_convolution(x,b,h)
    end if

end function wglTMDPDF_OPE_tw2

!!! this is convolution with twist3 PDF
function wglTMDPDF_OPE_tw3(x,b,h)
    real(dp),dimension(-5:5)::wglTMDPDF_OPE_tw3
    real(dp),intent(in)::x,b
    integer,intent(in)::h

    !!!! test for boundaries is done in wglTMDPDF_lowScale5 (on the entry to this procedure)

    !!! computation
    if(useGridTW3) then
        wglTMDPDF_OPE_tw3=mainGridTw3%Extract(x,b,h)
    else
        wglTMDPDF_OPE_tw3=tw3_convolution(x,b,h)
    end if

end function wglTMDPDF_OPE_tw3

!!!!!!!!!! ------------------------ SUPPORTING ROUTINES --------------------------------------
!!! This subroutine forces reconstruction of the grid for tw2 part (if gridding is ON)
subroutine wglTMDPDF_OPE_tw2_resetGrid()
    if(useGridTW2) then
        if(outputLevel>1) write(*,*) 'arTeMiDe ',moduleName,':  Grid Reset for tw2. with c4=',c4_tw2_global
        call mainGridTw2%MakeGrid(tw2_convolution)!
    end if
end subroutine wglTMDPDF_OPE_tw2_resetGrid

!!! This subroutine forces reconstruction of the grid for tw3 part (if gridding is ON)
subroutine wglTMDPDF_OPE_tw3_resetGrid()
    !!!!! not implemented yet
end subroutine wglTMDPDF_OPE_tw3_resetGrid

!! call QCDinput to change the PDF replica number for tw2 part
!! unset the grid, since it should be recalculated for different PDF replica.
subroutine wglTMDPDF_OPE_tw2_SetPDFreplica(rep,hadron)
    integer,intent(in):: rep,hadron
    logical::newPDF

    call QCDinput_SethPDFreplica(rep,hadron,newPDF)
    if(newPDF) then
        call wglTMDPDF_OPE_tw2_resetGrid()
    else
        if(outputLevel>1) write(*,"('arTeMiDe ',A,':  replica of PDF (',I4,') is the same as the used one. Nothing is done!')") &
        moduleName, rep
    end if

end subroutine wglTMDPDF_OPE_tw2_SetPDFreplica

!! call QCDinput to change the PDF replica number for tw3 part
!! unset the grid, since it should be recalculated for different PDF replica.
subroutine wglTMDPDF_OPE_tw3_SetPDFreplica(rep,hadron)
    integer,intent(in):: rep,hadron
    !!!!! not implemented yet
end subroutine wglTMDPDF_OPE_tw3_SetPDFreplica

!!!! this routine sets the variations of scales for tw2 part
!!!! it is used for the estimation of errors
subroutine wglTMDPDF_OPE_tw2_SetScaleVariation(c4_in)
    real(dp),intent(in)::c4_in
    if(c4_in<0.1d0 .or. c4_in>10.d0) then
        if(outputLevel>0) write(*,*) WarningString('variation in c4 (tw2) is enourmous. c4 is set to 2',moduleName)
        c4_tw2_global=2d0
        call wglTMDPDF_OPE_tw2_resetGrid()
    else if(abs(c4_in-c4_tw2_global)<toleranceGEN) then
        if(outputLevel>1) write(*,*) color('wglTMDPDF: c4-variation is ignored (tw2). c4='//numToStr(c4_tw2_global),c_yellow)
    else
        c4_tw2_global=c4_in
        if(outputLevel>1) write(*,*) color('wglTMDPDF: set scale variations c4 (tw2) as:'//numToStr(c4_tw2_global),c_yellow)
        call wglTMDPDF_OPE_tw2_resetGrid()
    end if
end subroutine wglTMDPDF_OPE_tw2_SetScaleVariation

!!!! this routine sets the variations of scales  for tw3 part
!!!! it is used for the estimation of errors
subroutine wglTMDPDF_OPE_tw3_SetScaleVariation(c4_in)
    real(dp),intent(in)::c4_in
    if(c4_in<0.1d0 .or. c4_in>10.d0) then
        if(outputLevel>0) write(*,*) WarningString('variation in c4 (tw3) is enourmous. c4 is set to 2',moduleName)
        c4_tw3_global=2d0
        call wglTMDPDF_OPE_tw3_resetGrid()
    else if(abs(c4_in-c4_tw3_global)<toleranceGEN) then
        if(outputLevel>1) write(*,*) color('wglTMDPDF: c4-variation is ignored (tw3). c4='//numToStr(c4_tw3_global),c_yellow)
    else
        c4_tw3_global=c4_in
        if(outputLevel>1) write(*,*) color('wglTMDPDF: set scale variations c4 (tw3) as:'//numToStr(c4_tw3_global),c_yellow)
        call wglTMDPDF_OPE_tw3_resetGrid()
    end if
end subroutine wglTMDPDF_OPE_tw3_SetScaleVariation

end module wglTMDPDF_OPE
