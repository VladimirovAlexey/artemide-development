!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	Evaluation of the liniearly polarized gluon TMD PDF at low normalization point in zeta-prescription.
!	
!	if you use this module please, quote 1706.01473
!
!	20.08.2023  Implementation in ver.3.0
!
!				A.Vladimirov (20.08.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lpTMDPDF
use aTMDe_Numerics
use IO_functions
use QCDinput
use TMDR
use lpTMDPDF_OPE
use lpTMDPDF_model

implicit none
!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v3.00"
character (len=8),parameter :: moduleName="lpTMDPDF"
!Last appropriate version of constants-file
integer,parameter::inputver=30

!--------------------------------Working variables-----------------------------------------------
!--- general
logical:: started=.false.
!! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
integer::outputLevel=2
!! variable that count number of WRNING mesagges. In order not to spam too much
integer::messageTrigger=6

!!! the length and array of NP parameters
integer::lambdaNPlength
real(dp),dimension(:),allocatable::lambdaNP
real(dp)::BMAX_ABS=100._dp !!! for large values of b returns 0
real(dp)::toleranceGEN !!! tolerance general

integer :: messageCounter

!!! General parameters
logical::includeGluon=.true.    !! gluons included/non-included (TRUE for lp TMDPDF!)
integer::numOfHadrons=1         !! total number of hadrons to compute

!!-----------------------------------------------Public interface---------------------------------------------------

public::lpTMDPDF_Initialize,lpTMDPDF_IsInitialized,lpTMDPDF_SetScaleVariation,lpTMDPDF_SetPDFreplica
public::lpTMDPDF_SetLambdaNP,lpTMDPDF_CurrentLambdaNP
public::lpTMDPDF_lowScale5

interface lpTMDPDF_inB
    module procedure lpTMDPDF_opt,lpTMDPDF_ev
end interface


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Interface subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function lpTMDPDF_IsInitialized()
    logical::lpTMDPDF_IsInitialized
    lpTMDPDF_IsInitialized=started
end function lpTMDPDF_IsInitialized

!! Initialization of the package
subroutine lpTMDPDF_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path
    logical::initRequired
    integer::FILEver

    if(started) return

    if(present(prefix)) then
        path=trim(adjustl(prefix))//trim(adjustr(file))
    else
        path=trim(adjustr(file))
    end if

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
        CLOSE (51, STATUS='KEEP')
        stop
    end if

    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>1) write(*,*) '--------------------------------------------- '
    if(outputLevel>1) write(*,*) 'artemide.',moduleName,version,': initialization started ... '

    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger

    call MoveTO(51,'*11  ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>1) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        CLOSE (51, STATUS='KEEP')
        return
    end if

!-------------general parameters
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) includeGluon

    call MoveTO(51,'*p2  ')
    read(51,*) numOfHadrons

    !-------------parameters of NP model
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) lambdaNPlength

    call MoveTO(51,'*p2  ')
    read(51,*) BMAX_ABS


    if(lambdaNPlength<=0) then
        write(*,*) ErrorString(&
        'Initialize: number of non-pertrubative parameters should be >=1. Check the constants-file. Evaluation STOP',moduleName)
            CLOSE (51, STATUS='KEEP')
        stop
    end if

    !!!!! ---- parameters of numerical evaluation
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p2  ')
    read(51,*) toleranceGEN

    CLOSE (51, STATUS='KEEP')

    if(.not.includeGluon) then
        write(*,*) color('INCONSITENCY: lpTMDPDF should include gluon',c_red)
    end if

    if(outputLevel>2 .and. includeGluon) write(*,'(A)') ' ... gluons are included'
    if(outputLevel>2 .and. .not.includeGluon) write(*,'(A)') ' ... gluons are not included'
    if(outputLevel>2) write(*,'(A,I3)') ' Number of hadrons to be considered =',numOfHadrons
    if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',lambdaNPlength
    if(outputLevel>2) write(*,'(A,F12.2)') ' Absolute maximum b      =',BMAX_ABS

    allocate(lambdaNP(1:lambdaNPlength))


    if(.not.lpTMDPDF_OPE_IsInitialized()) then
        if(outputLevel>2) write(*,*) '.. initializing lpTMDPDF_OPE (from ',moduleName,')'
        if(present(prefix)) then
            call lpTMDPDF_OPE_Initialize(file,prefix)
        else
            call lpTMDPDF_OPE_Initialize(file)
        end if
    end if

    call ModelInitialization(lambdaNPlength)
    if(outputLevel>0) write(*,*) color('----- arTeMiDe.lpTMDPDF_model : .... initialized',c_green)

    started=.true.
    messageCounter=0

    if(outputLevel>0) write(*,*) color('----- arTeMiDe.lpTMDPDF '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '
end subroutine lpTMDPDF_Initialize

!!!!!!!!!! ------------------------ SUPPORINTG ROUTINES --------------------------------------
!!! update PDF replica
subroutine lpTMDPDF_SetPDFreplica(rep,hadron)
    integer,intent(in):: rep,hadron

    call lpTMDPDF_OPE_SetPDFreplica(rep,hadron)
end subroutine lpTMDPDF_SetPDFreplica

!!!! this routine set the variations of scales
subroutine lpTMDPDF_SetScaleVariation(c4_in)
    real(dp),intent(in)::c4_in
    call lpTMDPDF_OPE_SetScaleVariation(c4_in)
end subroutine lpTMDPDF_SetScaleVariation

!!!Sets the non-pertrubative parameters lambda
!!! carries additionl option to build the grid
!!! if need to build grid, specify the gluon required directive.
subroutine lpTMDPDF_SetLambdaNP(lambdaIN)
    real(dp),intent(in)::lambdaIN(:)
    integer::ll
    messageCounter=0

    ll=size(lambdaIN)
    if(ll/=lambdaNPlength) then
        if(outputLevel>0) write(*,"(A,I3,A,I3,')')")&
                WarningString('SetLambdaNP:length of lambdaNP(',moduleName),&
            ll,color(') does not match the required (',c_red),lambdaNPlength
        if(outputLevel>0) write(*,*)color('                NOTHING IS DONE!',c_red)
        return
    end if

    lambdaNP=lambdaIN
    call ModelUpdate(lambdaNP)

    if(outputLevel>2) write(*,*) 'arTeMiDe.',moduleName,': NPparameters reset = (',lambdaNP,')'
    call ModelUpdate(lambdaNP)
end subroutine lpTMDPDF_SetLambdaNP

!!! returns current value of NP parameters
function lpTMDPDF_CurrentLambdaNP()
    real(dp),dimension(1:lambdaNPlength)::lpTMDPDF_CurrentLambdaNP
    lpTMDPDF_CurrentLambdaNP=lambdaNP
end function lpTMDPDF_CurrentLambdaNP

!!!!!!!--------------------------- DEFINING ROUTINES ------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!======TO REMOVE
!!!!!!! the function that actually returns the lpTMDPDF!
function lpTMDPDF_lowScale5(x,bT,hadron)
  real(dp),dimension(-5:5)::lpTMDPDF_lowScale5
  real(dp),intent(in) :: x, bT
  integer,intent(in)::hadron

  !!! test boundaries
    if(x>1d0) then
        call Warning_Raise('Called x>1 (return 0). x='//numToStr(x),messageCounter,messageTrigger,moduleName)
        lpTMDPDF_lowScale5=0._dp
        return
    else if(x==1.d0) then !!! funny but sometimes FORTRAN can compare real numbers exactly
        lpTMDPDF_lowScale5=0._dp
        return
    else if(bT>BMAX_ABS) then
        lpTMDPDF_lowScale5=0._dp
        return
    else if(x<1d-12) then
        write(*,*) ErrorString('Called x<0. x='//numToStr(x)//' . Evaluation STOP',moduleName)
        stop
    else if(bT<0d0) then
        write(*,*) ErrorString('Called b<0. b='//numToStr(bT)//' . Evaluation STOP',moduleName)
        stop
    end if

    lpTMDPDF_lowScale5=lpTMDPDF_OPE_convolution(x,bT,hadron)*FNP(x,bT,hadron,lambdaNP)

    if(hadron<0) lpTMDPDF_lowScale5=lpTMDPDF_lowScale5(5:-5:-1)  !! it makes no sense, because it is only gluon... but still

end function lpTMDPDF_lowScale5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!======TO REMOVE

!!!!!!! the function that actually returns the lpTMDPDF optimal value
function lpTMDPDF_opt(x,bT,hadron)
  real(dp),dimension(-5:5)::lpTMDPDF_opt
  real(dp),intent(in) :: x, bT
  integer,intent(in)::hadron

  !!! test boundaries
    if(x>1d0) then
        call Warning_Raise('Called x>1 (return 0). x='//numToStr(x),messageCounter,messageTrigger,moduleName)
        lpTMDPDF_opt=0._dp
        return
     else if(x==1.d0) then !!! funny but sometimes FORTRAN can compare real numbers exactly
        lpTMDPDF_opt=0._dp
        return
    else if(bT>BMAX_ABS) then
        lpTMDPDF_opt=0._dp
        return
    else if(x<1d-12) then
        write(*,*) ErrorString('Called x<0. x='//numToStr(x)//' . Evaluation STOP',moduleName)
        stop
    else if(bT<0d0) then
        write(*,*) ErrorString('Called b<0. b='//numToStr(bT)//' . Evaluation STOP',moduleName)
        stop
    end if

    lpTMDPDF_opt=lpTMDPDF_OPE_convolution(x,bT,abs(hadron))*FNP(x,bT,abs(hadron),lambdaNP)

    if(hadron<0) lpTMDPDF_opt=lpTMDPDF_opt(5:-5:-1)

end function lpTMDPDF_opt
!
!!!!!!!! the function that actually returns the lpTMDPDF evolved to (mu,zeta) value
function lpTMDPDF_Ev(x,bt,muf,zetaf,hadron)
    real(dp)::lpTMDPDF_Ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: RkernelG

    if(includeGluon) then
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)

        lpTMDPDF_Ev=RkernelG*lpTMDPDF_opt(x,bT,hadron)

    else
        lpTMDPDF_Ev=0._dp
    end if


    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
    lpTMDPDF_Ev(5)=0_dp
    lpTMDPDF_Ev(-5)=0_dp
    end if
    if(muf<mCHARM) then
    lpTMDPDF_Ev(4)=0_dp
    lpTMDPDF_Ev(-4)=0_dp
    end if
end function lpTMDPDF_Ev

end module lpTMDPDF
