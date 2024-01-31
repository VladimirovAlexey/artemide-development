!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	Evaluation of the BoerMulders TMD PDF at low normalization point in zeta-prescription.
!
!	29.01.2024  Implementation in ver.3.0
!
!				A.Vladimirov (29.01.2024)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module BoerMuldersTMDPDF
use aTMDe_Numerics
use IO_functions
use BoerMuldersTMDPDF_OPE
use BoerMuldersTMDPDF_model

implicit none
!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v3.00"
character (len=17),parameter :: moduleName="BoerMuldersTMDPDF"
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

!!-----------------------------------------------Public interface---------------------------------------------------

public::BoerMuldersTMDPDF_Initialize,BoerMuldersTMDPDF_IsInitialized
public::BoerMuldersTMDPDF_SetScaleVariation_tw3
public::BoerMuldersTMDPDF_SetPDFreplica_tw3
public::BoerMuldersTMDPDF_SetLambdaNP,BoerMuldersTMDPDF_CurrentLambdaNP
public::BoerMuldersTMDPDF_lowScale5

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Interface subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function BoerMuldersTMDPDF_IsInitialized()
    logical::BoerMuldersTMDPDF_IsInitialized
    BoerMuldersTMDPDF_IsInitialized=started
end function BoerMuldersTMDPDF_IsInitialized

!! Initialization of the package
subroutine BoerMuldersTMDPDF_Initialize(file,prefix)
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

    call MoveTO(51,'*14   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>1) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        CLOSE (51, STATUS='KEEP')
        return
    end if

    !-------------parameters of NP model
    call MoveTO(51,'*B   ')
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

    if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',lambdaNPlength
    if(outputLevel>2) write(*,'(A,F12.2)') ' Absolute maximum b      =',BMAX_ABS

    allocate(lambdaNP(1:lambdaNPlength))

    !!!!! ---- parameters of numerical evaluation
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p2  ')
    read(51,*) toleranceGEN

    CLOSE (51, STATUS='KEEP') 

    if(.not.BoerMuldersTMDPDF_OPE_IsInitialized()) then
        if(outputLevel>2) write(*,*) '.. initializing BoerMuldersTMDPDF_OPE (from ',moduleName,')'
        if(present(prefix)) then
            call BoerMuldersTMDPDF_OPE_Initialize(file,prefix)
        else
            call BoerMuldersTMDPDF_OPE_Initialize(file)
        end if
    end if

    call ModelInitialization(lambdaNPlength)
    if(outputLevel>0) write(*,*) color('----- arTeMiDe.BoerMuldersTMDPDF_model : .... initialized',c_green)

    started=.true.
    messageCounter=0

    if(outputLevel>0) write(*,*) color('----- arTeMiDe.BoerMuldersTMDPDF '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '
end subroutine BoerMuldersTMDPDF_Initialize

!!!!!!!!!! ------------------------ SUPPORINTG ROUTINES --------------------------------------
!!! update PDF replica
subroutine BoerMuldersTMDPDF_SetPDFreplica_tw3(rep,hadron)
    integer,intent(in):: rep,hadron

    call BoerMuldersTMDPDF_OPE_tw3_SetPDFreplica(rep,hadron)
end subroutine BoerMuldersTMDPDF_SetPDFreplica_tw3

!!!! this routine set the variations of scales
subroutine BoerMuldersTMDPDF_SetScaleVariation_tw3(c4_in)
    real(dp),intent(in)::c4_in
    call BoerMuldersTMDPDF_OPE_tw3_SetScaleVariation(c4_in)
end subroutine BoerMuldersTMDPDF_SetScaleVariation_tw3

!!!Sets the non-pertrubative parameters lambda
!!! carries additionl option to build the grid
!!! if need to build grid, specify the gluon required directive.
subroutine BoerMuldersTMDPDF_SetLambdaNP(lambdaIN)
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
end subroutine BoerMuldersTMDPDF_SetLambdaNP

!!! returns current value of NP parameters
function BoerMuldersTMDPDF_CurrentLambdaNP()
    real(dp),dimension(1:lambdaNPlength)::BoerMuldersTMDPDF_CurrentLambdaNP
    BoerMuldersTMDPDF_CurrentLambdaNP=lambdaNP
end function BoerMuldersTMDPDF_CurrentLambdaNP

!!!!!!!--------------------------- DEFINING ROUTINES ------------------------------------------

!!!!!!! the function that actually returns the BoerMuldersTMDPDF!
function BoerMuldersTMDPDF_lowScale5(x,bT,hadron)
  real(dp),dimension(-5:5)::BoerMuldersTMDPDF_lowScale5
  real(dp),intent(in) :: x, bT
  integer,intent(in)::hadron

  !!! test boundaries
    if(x>1d0) then
        call Warning_Raise('Called x>1 (return 0). x='//numToStr(x),messageCounter,messageTrigger,moduleName)
        BoerMuldersTMDPDF_lowScale5=0._dp
        return
    else if(x==1.d0) then !!! funny but sometimes FORTRAN can compare real numbers exactly
        BoerMuldersTMDPDF_lowScale5=0._dp
        return
    else if(bT>BMAX_ABS) then
        BoerMuldersTMDPDF_lowScale5=0._dp
        return
    else if(x<1d-12) then
        write(*,*) ErrorString('Called x<0. x='//numToStr(x)//' . Evaluation STOP',moduleName)
        stop
    else if(bT<0d0) then
        write(*,*) ErrorString('Called b<0. b='//numToStr(bT)//' . Evaluation STOP',moduleName)
        stop
    end if

    BoerMuldersTMDPDF_lowScale5=BoerMuldersTMDPDF_OPE_tw3_convolution(x,bT,abs(hadron))*FNP(x,bT,abs(hadron),lambdaNP)

    if(hadron<0) BoerMuldersTMDPDF_lowScale5=BoerMuldersTMDPDF_lowScale5(5:-5:-1)

end function BoerMuldersTMDPDF_lowScale5

end module BoerMuldersTMDPDF