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
use aTMDe_IO
use aTMDe_Ogata
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
type(Warning_OBJ)::Warning_Handler

!!! the length and array of NP parameters
integer::lambdaNPlength
real(dp),dimension(:),allocatable::lambdaNP
real(dp)::BMAX_ABS=100._dp !!! for large values of b returns 0
real(dp)::toleranceGEN !!! tolerance general

!!!------------------------------ General parameters----------------------------------------------
logical::includeGluon=.true.    !! gluons included/non-included (TRUE for lp TMDPDF!)
integer::numOfHadrons=1         !! total number of hadrons to compute
real(dp)::TMDmass=1._dp         !! mass parameter used as mass-scale

!!!------------------------------ Parameters of transform to KT-space -------------------------------------------

integer,parameter::TMDtypeN=2 !!!!! this is the order of Bessel-transform (IT IS STRICT FOR TMD)
real(dp)::kT_FREEZE=0.0001_dp  !!!!! parameter of freezing the low-kT-value

type(OgataIntegrator)::Hankel

!!!------------------------------ Parameters of transform to TMM -------------------------------------------

real(dp)::muTMM_min=0.8_dp  !!!!! minimal mu

!!-----------------------------------------------Public interface---------------------------------------------------

public::lpTMDPDF_Initialize,lpTMDPDF_IsInitialized,lpTMDPDF_SetScaleVariation,lpTMDPDF_SetPDFreplica
public::lpTMDPDF_SetLambdaNP,lpTMDPDF_CurrentLambdaNP
public::lpTMDPDF_inB,lpTMDPDF_inKT,lpTMDPDF_TMM_G,lpTMDPDF_TMM_X

interface lpTMDPDF_inB
    module procedure TMD_opt,TMD_ev
end interface

interface lpTMDPDF_inKT
    module procedure Fourier_opt,Fourier_ev
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
integer::FILEver,messageTrigger
real(dp)::hOGATA,toleranceOGATA
real(dp)::hOGATA_TMM,toleranceOGATA_TMM

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
    ERROR STOP
end if

call MoveTO(51,'*p2  ')
read(51,*) outputLevel
if(outputLevel>1) write(*,*) '--------------------------------------------- '
if(outputLevel>1) write(*,*) 'artemide.',moduleName,version,': initialization started ... '

call MoveTO(51,'*p3  ')
read(51,*) messageTrigger

call MoveTO(51,'*B   ')
call MoveTO(51,'*p2  ')
read(51,*) TMDmass

    !! TMDR
call MoveTO(51,'*3   ')
call MoveTO(51,'*p1  ')
read(51,*) initRequired
if(.not.initRequired) then
    write(*,*) ErrorString('TMDR module MUST be included.',moduleName)
    write(*,*) ErrorString('Check initialization-file. Evaluation stop.',moduleName)
    CLOSE (51, STATUS='KEEP')
    ERROR STOP
end if

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
    'Initialize: number of non-perturbative parameters should be >=1. Check the constants-file. Evaluation STOP',moduleName)
        CLOSE (51, STATUS='KEEP')
    ERROR STOP
end if

!!!!! ---- parameters of numerical evaluation
call MoveTO(51,'*D   ')
call MoveTO(51,'*p2  ')
read(51,*) toleranceGEN


!!!!! ---- parameters of KT-transformation
call MoveTO(51,'*F   ')
call MoveTO(51,'*p1  ')
read(51,*) toleranceOGATA
call MoveTO(51,'*p2  ')
read(51,*) hOGATA
call MoveTO(51,'*p3  ')
read(51,*) kT_FREEZE

!!!!! ---- parameters of TMM-transformation
call MoveTO(51,'*G   ')
call MoveTO(51,'*p1  ')
read(51,*) toleranceOGATA_TMM
call MoveTO(51,'*p2  ')
read(51,*) hOGATA_TMM
call MoveTO(51,'*p3  ')
read(51,*) muTMM_min

CLOSE (51, STATUS='KEEP')
Warning_Handler=Warning_OBJ(moduleName=moduleName,messageCounter=0,messageTrigger=messageTrigger)

if(.not.includeGluon) then
    write(*,*) color('INCONSITENCY: lpTMDPDF should include gluon',c_red)
end if

if(outputLevel>2 .and. includeGluon) write(*,'(A)') ' ... gluons are included'
if(outputLevel>2 .and. .not.includeGluon) write(*,'(A)') ' ... gluons are not included'
if(outputLevel>2) write(*,'(A,I3)') ' Number of hadrons to be considered =',numOfHadrons
if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',lambdaNPlength
if(outputLevel>2) write(*,'(A,F12.2)') ' Absolute maximum b      =',BMAX_ABS

allocate(lambdaNP(1:lambdaNPlength))

Hankel=OgataIntegrator(moduleName,outputLevel,TMDtypeN, toleranceOGATA,hOGATA,TMDmass)

if(.not.TMDR_IsInitialized()) then
    if(outputLevel>2) write(*,*) '.. initializing TMDR (from ',moduleName,')'
    if(present(prefix)) then
        call TMDR_Initialize(file,prefix)
    else
        call TMDR_Initialize(file)
    end if
end if

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
    call Warning_Handler%Reset()

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

!!!!! the names are neutral because these procedures are feed to Fourier transform. And others universal sub programs.

!!!!!!! the function that actually returns the lpTMDPDF optimal value
function TMD_opt(x,bT,hadron)
  real(dp),dimension(-5:5)::TMD_opt
  real(dp),intent(in) :: x, bT
  integer,intent(in)::hadron

  !!! test boundaries
    if(x>1d0) then
        call Warning_Handler%WarningRaise('Called x>1 (return 0). x='//numToStr(x))
        TMD_opt=0._dp
        return
     else if(x==1.d0) then !!! funny but sometimes FORTRAN can compare real numbers exactly
        TMD_opt=0._dp
        return
    else if(bT>BMAX_ABS) then
        TMD_opt=0._dp
        return
    else if(x<1d-12) then
        ERROR STOP ErrorString('Called x<0. x='//numToStr(x)//' . Evaluation STOP',moduleName)
        stop
    else if(bT<0d0) then
        ERROR STOP ErrorString('Called b<0. b='//numToStr(bT)//' . Evaluation STOP',moduleName)
        stop
    end if

    TMD_opt=lpTMDPDF_OPE_convolution(x,bT,abs(hadron))*FNP(x,bT,abs(hadron),lambdaNP)

    if(hadron<0) TMD_opt=TMD_opt(5:-5:-1)

end function TMD_opt
!
!!!!!!!! the function that actually returns the lpTMDPDF evolved to (mu,zeta) value
function TMD_ev(x,bt,muf,zetaf,hadron)
    real(dp)::TMD_ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: RkernelG

    if(includeGluon) then
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)

        TMD_ev=RkernelG*TMD_opt(x,bT,hadron)

    else
        TMD_ev=0._dp
    end if


    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
    TMD_ev(5)=0_dp
    TMD_ev(-5)=0_dp
    end if
    if(muf<mCHARM) then
    TMD_ev(4)=0_dp
    TMD_ev(-4)=0_dp
    end if
end function TMD_ev

!!!!!!!! TMM G_{n,n} at (x,mu)
function lpTMDPDF_TMM_G(x,mu,hadron)
    real(dp)::lpTMDPDF_TMM_G(-5:5)
    real(dp),intent(in):: x,mu
    integer,intent(in)::hadron

    if(mu<muTMM_min) then
        write(*,*) ErrorString("ERROR in KT-moment computation. mu<mu_min:",moduleName),muTMM_min
        error stop
    end if

    lpTMDPDF_TMM_G=Hankel%Moment_G(F,mu)

    !!! forcefully set =0 below threshold
    if(mu<mBOTTOM) then
    lpTMDPDF_TMM_G(5)=0_dp
    lpTMDPDF_TMM_G(-5)=0_dp
    end if
    if(mu<mCHARM) then
    lpTMDPDF_TMM_G(4)=0_dp
    lpTMDPDF_TMM_G(-4)=0_dp
    end if
contains
    function F(b)
    real(dp),dimension(-5:5)::F
    real(dp),intent(in)::b
    F=TMD_opt(x,b,hadron)
    end function F
end function lpTMDPDF_TMM_G

!!!!!!!! TMM G_{n+1,n} at (x,mu)
function lpTMDPDF_TMM_X(x,mu,hadron)
    real(dp)::lpTMDPDF_TMM_X(-5:5)
    real(dp),intent(in):: x,mu
    integer,intent(in)::hadron

    if(mu<muTMM_min) then
        write(*,*) ErrorString("ERROR in KT-moment computation. mu<mu_min:",moduleName),muTMM_min
        error stop
    end if

    lpTMDPDF_TMM_X=Hankel%Moment_X(F,mu)

    !!! forcefully set =0 below threshold
    if(mu<mBOTTOM) then
    lpTMDPDF_TMM_X(5)=0_dp
    lpTMDPDF_TMM_X(-5)=0_dp
    end if
    if(mu<mCHARM) then
    lpTMDPDF_TMM_X(4)=0_dp
    lpTMDPDF_TMM_X(-4)=0_dp
    end if
contains
    function F(b)
    real(dp),dimension(-5:5)::F
    real(dp),intent(in)::b
    F=TMD_opt(x,b,hadron)
    end function F
end function lpTMDPDF_TMM_X

!!!--------------------------------------------------------------------------------------
!!!------------------------------------------FOURIER-------------------------------------
!!!--------------------------------------------------------------------------------------
!!! It evaluates the integral for the transformation to the kT-space
!!! int_0^infty   b db/2pi  J_num(b qT) F1  (b/qT)^num M^{2num}/num!
!!! the integration is made by the class aTMDe_Ogata
function Fourier_ev(x,qT_in,mu,zeta,hadron)
real(dp),intent(in)::x,mu,zeta,qT_in
integer,intent(in)::hadron
real(dp)::Fourier_ev(-5:5)

real(dp)::qT
if(qT_in<kT_FREEZE) then
    qT=kT_FREEZE
else
    qT=qT_in
end if

Fourier_ev=Hankel%TransformTMD(F,qT)

contains
function F(b)
real(dp),dimension(-5:5)::F
real(dp),intent(in)::b
F=TMD_ev(x,b,mu,zeta,hadron)
end function F

end function Fourier_ev

!!! It evaluates the integral for the transformation to the kT-space
!!! int_0^infty   b db/2pi  J_num(b qT) F1  (b/qT)^num M^{2num}/num!
!!! the integration is made by the class aTMDe_Ogata
function Fourier_opt(x,qT_in,hadron)
real(dp),intent(in)::x,qT_in
integer,intent(in)::hadron
real(dp)::Fourier_opt(-5:5)

real(dp)::qT
if(qT_in<kT_FREEZE) then
    qT=kT_FREEZE
else
    qT=qT_in
end if

Fourier_opt=Hankel%TransformTMD(F,qT)

contains
function F(b)
real(dp),dimension(-5:5)::F
real(dp),intent(in)::b
F=TMD_opt(x,b,hadron)
end function F

end function Fourier_opt

end module lpTMDPDF
