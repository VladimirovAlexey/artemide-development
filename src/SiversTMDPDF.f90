!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.05
!
!	Evaluation of the Sivers TMD PDF at low normalization point in zeta-prescription.
!
!	21.08.2023  Implementation in ver.3.0
!	23.06.2026  Removed old INCLUDE-dependencies + updated of KT-part
!   06.07.2026  Changing to a new h-enumeration scheme
!
!				A.Vladimirov (21.08.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module SiversTMDPDF
use aTMDe_Numerics
use aTMDe_IO
use aTMDe_Ogata
use aTMDe_Levin
use aTMDe_ktGrid
use aTMDe_hDef
use QCDinput
use TMDR
use SiversTMDPDF_OPE
use SiversTMDPDF_model

implicit none
!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v3.05"
character (len=12),parameter :: moduleName="SiversTMDPDF"
!Last appropriate version of constants-file
integer,parameter::inputver=41

!--------------------------------Working variables-----------------------------------------------
!--- general
logical:: started=.false.
!! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
!! 3=FULL output
integer::outputLevel=2
type(Warning_OBJ)::Warning_Handler

!!! the length and array of NP parameters
integer::lambdaNPlength
real(dp),dimension(:),allocatable::lambdaNP
real(dp)::BMAX_ABS=100._dp !!! for large values of b returns 0
real(dp)::toleranceGEN !!! tolerance general

!!!------------------------------ General parameters----------------------------------------------
logical::includeGluon=.false.   !! gluons included/excluded
integer::numOfHadrons=1         !! total number of hadrons to compute
real(dp)::TMDmass=1._dp         !! mass parameter used as a mass-scale

!!!------------------------------ Parameters of transform to KT-space -------------------------------------------

integer,parameter::TMDtypeN=1 !!!!! this is the order of Bessel-transform (IT IS STRICT FOR TMD)

type(OgataIntegrator)::HankelbyOGATA
real(dp),parameter::kT_FREEZE=0.0001_dp  !!!!! parameter of freezing the low-kT-value in OGATA (it is not dynamic in the current implementation)

!!!------------------------------ Parameters of transform in KT space and KT-grid ------------------------------
logical::makeGrid_inKT
type(ktGrid)::mainKTGrid
type(LevinIntegrator)::HankelbyLEVIN
!!!------------------------------ Parameters of transform to TMM -------------------------------------------

real(dp)::muTMM_min=0.8_dp  !!!!! minimal mu

!!-----------------------------------------------Public interface---------------------------------------------------

public::SiversTMDPDF_Initialize,SiversTMDPDF_IsInitialized
public::SiversTMDPDF_SetScaleVariation_tw3
public::SiversTMDPDF_SetPDFreplica_tw3
public::SiversTMDPDF_SetLambdaNP,SiversTMDPDF_CurrentLambdaNP
public::SiversTMDPDF_inB,SiversTMDPDF_inKT,SiversTMDPDF_TMM_G,SiversTMDPDF_TMM_X

interface SiversTMDPDF_inB
    module procedure TMD_opt,TMD_ev
end interface

interface SiversTMDPDF_inKT
    module procedure TMD_opt_inKT,TMD_ev_inKT, TMD_grid_inKT
end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Interface subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function SiversTMDPDF_IsInitialized()
    logical::SiversTMDPDF_IsInitialized
    SiversTMDPDF_IsInitialized=started
end function SiversTMDPDF_IsInitialized

!! Initialization of the package
subroutine SiversTMDPDF_Initialize(file,prefix)
character(len=*)::file
character(len=*),optional::prefix
character(len=:),allocatable::path
logical::initRequired
integer::FILEver,messageTrigger
real(dp)::hOGATA_TMM,toleranceOGATA_TMM

if(started) return

if(present(prefix)) then
    path=trim(adjustl(prefix))//trim(adjustl(file))
else
    path=trim(adjustl(file))
end if

OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")

call MoveTO(51,'*0   ')
call MoveTO(51,'*A   ')
call MoveTO(51,'*p1  ')
read(51,*) FILEver
if(FILEver<inputver) then
    write(*,*) 'artemide.'//trim(moduleName)//': const-file version is too old.'
    write(*,*) '		     Update the const-file with artemide.setup'
    write(*,*) '  '
    CLOSE (51, STATUS='KEEP')
    error stop
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
    error stop
end if

call MoveTO(51,'*12  ')
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
error stop
end if

!!!!! ---- parameters of numerical evaluation
call MoveTO(51,'*D   ')
call MoveTO(51,'*p2  ')
read(51,*) toleranceGEN

!!!!! ---- parameters of kT-transform and grid
call MoveTO(51,'*F   ')
call MoveTO(51,'*p1  ')
read(51,*) makeGrid_inKT

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

if(outputLevel>2 .and. includeGluon) write(*,'(A)') ' ... gluons are included'
if(outputLevel>2 .and. .not.includeGluon) write(*,'(A)') ' ... gluons are not included'
if(outputLevel>2) write(*,'(A,I3)') ' Number of hadrons to be considered =',numOfHadrons
if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',lambdaNPlength
if(outputLevel>2) write(*,'(A,F12.2)') ' Absolute maximum b      =',BMAX_ABS

allocate(lambdaNP(1:lambdaNPlength))

HankelbyLEVIN=LevinIntegrator(path,'*12  ','*F   ',moduleName,outputLevel,TMDtypeN)
if(makeGrid_inKT) then
        mainKTGrid=ktGrid(path,'*12   ','*F   ',numOfHadrons,includeGluon,moduleName,outputLevel)
end if

!!!!!! TODO: fix the minimal value of KT
HankelbyOGATA=OgataIntegrator(moduleName,outputLevel,TMDtypeN, toleranceOGATA_TMM,hOGATA_TMM,TMDmass,kT_FREEZE)

if(.not.TMDR_IsInitialized()) then
    if(outputLevel>2) write(*,*) '.. initializing TMDR (from ',moduleName,')'
    if(present(prefix)) then
        call TMDR_Initialize(file,prefix)
    else
        call TMDR_Initialize(file)
    end if
end if

if(.not.SiversTMDPDF_OPE_IsInitialized()) then
    if(outputLevel>2) write(*,*) '.. initializing SiversTMDPDF_OPE (from ',moduleName,')'
    if(present(prefix)) then
        call SiversTMDPDF_OPE_Initialize(file,prefix)
    else
        call SiversTMDPDF_OPE_Initialize(file)
    end if
end if

call ModelInitialization(lambdaNPlength)
if(outputLevel>0) write(*,*) color('----- arTeMiDe.SiversTMDPDF_model : .... initialized',c_green)

started=.true.

if(outputLevel>0) write(*,*) color('----- arTeMiDe.SiversTMDPDF '//trim(version)//': .... initialized',c_green)
if(outputLevel>1) write(*,*) ' '
end subroutine SiversTMDPDF_Initialize

!!!!!!!!!! ------------------------ SUPPORTING ROUTINES --------------------------------------
!!! update PDF replica
subroutine SiversTMDPDF_SetPDFreplica_tw3(rep,hadron)
    integer,intent(in):: rep,hadron

    call SiversTMDPDF_OPE_tw3_SetPDFreplica(rep,hadron)
end subroutine SiversTMDPDF_SetPDFreplica_tw3

!!!! this routine sets the variations of scales
subroutine SiversTMDPDF_SetScaleVariation_tw3(c4_in)
    real(dp),intent(in)::c4_in
    call SiversTMDPDF_OPE_tw3_SetScaleVariation(c4_in)
end subroutine SiversTMDPDF_SetScaleVariation_tw3

!!!Sets the non-perturbative parameters lambda
subroutine SiversTMDPDF_SetLambdaNP(lambdaIN)
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

    if(makeGrid_inKT) call updateGrid_inKT()
end subroutine SiversTMDPDF_SetLambdaNP

!!! returns current value of NP parameters
function SiversTMDPDF_CurrentLambdaNP()
    real(dp),dimension(1:lambdaNPlength)::SiversTMDPDF_CurrentLambdaNP
    SiversTMDPDF_CurrentLambdaNP=lambdaNP
end function SiversTMDPDF_CurrentLambdaNP

!!!!! Function that constructs the grid in KT-space
!!!!! I take the TMD and send it to Fourier, and then resend it to the grid...
subroutine updateGrid_inKT()
real(dp)::Q,x
integer::h

call mainKTGrid%MakeGrid(toGrid)

contains

function toFourier(b)
    real(dp),dimension(-5:5) :: toFourier
    real(dp), intent(in) ::b

    toFourier=TMD_ev(x,b,Q,Q**2,abs(h))
end function toFourier

function toGrid(x_in,Q_in,h_in,arraySize1,arraySize2)
    integer,intent(in)::arraySize1,arraySize2
    real(dp),dimension(1:arraySize1,0:arraySize2,-5:5) :: toGrid
    real(dp), intent(in) ::x_in,Q_in
    integer,intent(in)::h_in

    x=x_in
    Q=Q_in
    h=h_in

    toGrid=HankelbyLEVIN%Fourier_array(toFourier)
end function toGrid

end subroutine updateGrid_inKT

!!!!!!!--------------------------- DEFINING ROUTINES ------------------------------------------

!!!!! the names are neutral because these procedures are fed to Fourier transform and other universal subprograms.

!!!!!!! the function that actually returns the SiversTMDPDF!
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
    else if(x<=0) then
        error stop ErrorString('Called x<0. x='//numToStr(x)//' . Evaluation STOP',moduleName)
    else if(bT<0d0) then
        error stop ErrorString('Called b<0. b='//numToStr(bT)//' . Evaluation STOP',moduleName)
    end if

    if(abs(hadron)<10) then
        !!! This is usual hadron
        TMD_opt=SiversTMDPDF_OPE_tw3(x,bT,abs(hadron))*FNP(x,bT,abs(hadron),lambdaNP)
        if(hadron<0) TMD_opt=TMD_opt(5:-5:-1)
    else
        !!! This is an unusual hadron (all unusual hadrons are derived from hadron 1)
        TMD_opt=SiversTMDPDF_OPE_tw3(x,bT,1)*FNP(x,bT,1,lambdaNP)
        TMD_opt=PDFforH(TMD_opt,hadron)
    end if

end function TMD_opt

!!!!!!!! the function that actually returns the SiversTMDPDF evolved to (mu,zeta) value
function TMD_ev(x,bt,muf,zetaf,hadron)
    real(dp)::TMD_ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: Rkernel,RkernelG

    if(bT>BMAX_ABS) then
      TMD_ev=0._dp
      return
    end if

    if(includeGluon) then
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)

        TMD_ev=TMD_opt(x,bT,hadron)*&
            (/Rkernel,Rkernel,Rkernel,Rkernel,Rkernel,RkernelG,Rkernel,Rkernel,Rkernel,Rkernel,Rkernel/)

    else
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
        TMD_ev=Rkernel*TMD_opt(x,bT,hadron)
    end if


    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
    TMD_ev(5)=0._dp
    TMD_ev(-5)=0._dp
    end if
    if(muf<mCHARM) then
    TMD_ev(4)=0._dp
    TMD_ev(-4)=0._dp
    end if

end function TMD_ev

!!!!!!! the function that actually returns the optimal value
function TMD_opt_inKT(x,kT,hadron)
  real(dp),dimension(-5:5)::TMD_opt_inKT
  real(dp),intent(in) :: x, kT
  integer,intent(in)::hadron

  !!! test boundaries
    if(x>1d0) then
        call Warning_Handler%WarningRaise('Called x>1 (return 0). x='//numToStr(x))
        TMD_opt_inKT=0._dp
        return
     else if(x==1.d0) then !!! funny but sometimes FORTRAN can compare real numbers exactly
        TMD_opt_inKT=0._dp
        return
    else if(x<=0) then
        error stop ErrorString('Called x<0. x='//numToStr(x)//' . Evaluation STOP',moduleName)
    else if(kT<0d0) then
        error stop ErrorString('Called kT<0. kT='//numToStr(kT)//' . Evaluation STOP',moduleName)
    end if

    TMD_opt_inKT=HankelbyLEVIN%Fourier_atPoint(toFourier,kT)

    if(hadron<0) TMD_opt_inKT=TMD_opt_inKT(5:-5:-1)

    contains

    function toFourier(b)
        real(dp),dimension(-5:5) :: toFourier
        real(dp), intent(in) ::b
        toFourier=TMD_opt(x,b,abs(hadron))
    end function toFourier

end function TMD_opt_inKT
!
!!!!!!!! the function that actually returns the TMD evolved to (mu,zeta) value
function TMD_ev_inKT(x,kT,muf,zetaf,hadron)
    real(dp)::TMD_ev_inKT(-5:5)
    real(dp),intent(in):: x,kT,muf,zetaf
    integer,intent(in)::hadron

    if(x>=1._dp) then
      TMD_ev_inKT=0._dp
      return
    end if

    TMD_ev_inKT=HankelbyLEVIN%Fourier_atPoint(toFourier,kT)

    if(hadron<0) TMD_ev_inKT=TMD_ev_inKT(5:-5:-1)

    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
    TMD_ev_inKT(5)=0._dp
    TMD_ev_inKT(-5)=0._dp
    end if
    if(muf<mCHARM) then
    TMD_ev_inKT(4)=0._dp
    TMD_ev_inKT(-4)=0._dp
    end if

contains

function toFourier(b)
    real(dp),dimension(-5:5) :: toFourier
    real(dp), intent(in) ::b
    toFourier=TMD_ev(x,b,muf,zetaf,abs(hadron))
end function toFourier

end function TMD_ev_inKT


!!!!!!!! the function that actually returns the TMD evolved to (mu,mu^2) value
!!!!!!! This is exactly what is stored in the grid. So if the grid is built, extract from it  .
function TMD_grid_inKT(x,kT,mu,hadron)
    real(dp)::TMD_grid_inKT(-5:5)
    real(dp),intent(in):: x,kT,mu
    integer,intent(in)::hadron

    if(makeGrid_inKT) then
        if(abs(hadron)<10) then
            !!!! this is usual hadron
            TMD_grid_inKT=mainKTGrid%Extract(x,kT,mu,abs(hadron))
            if(hadron<0) TMD_grid_inKT=TMD_grid_inKT(5:-5:-1)
        else
            !!!! This is unusual hadron (all unusual hadrons are made from hadron 1)
            TMD_grid_inKT=mainKTGrid%Extract(x,kT,mu,1)
            TMD_grid_inKT=PDFforH(TMD_grid_inKT,hadron)
        end if
    else
        TMD_grid_inKT=TMD_ev_inKT(x,kT,mu,mu**2,hadron)
    end if

end function TMD_grid_inKT

!!!!!!!! TMM G_{n,n} at (x,mu)
function SiversTMDPDF_TMM_G(x,mu,hadron)
    real(dp)::SiversTMDPDF_TMM_G(-5:5)
    real(dp),intent(in):: x,mu
    integer,intent(in)::hadron

    if(mu<muTMM_min) then
        write(*,*) ErrorString("ERROR in KT-moment computation. mu<mu_min:",moduleName),muTMM_min
        error stop
    end if

    SiversTMDPDF_TMM_G=HankelbyOGATA%Moment_G(F,mu)

    !!! forcefully set =0 below threshold
    if(mu<mBOTTOM) then
    SiversTMDPDF_TMM_G(5)=0._dp
    SiversTMDPDF_TMM_G(-5)=0._dp
    end if
    if(mu<mCHARM) then
    SiversTMDPDF_TMM_G(4)=0._dp
    SiversTMDPDF_TMM_G(-4)=0._dp
    end if
contains
    function F(b)
    real(dp),dimension(-5:5)::F
    real(dp),intent(in)::b
    F=TMD_opt(x,b,hadron)
    end function F

end function SiversTMDPDF_TMM_G

!!!!!!!! TMM G_{n+1,n} at (x,mu)
function SiversTMDPDF_TMM_X(x,mu,hadron)
    real(dp)::SiversTMDPDF_TMM_X(-5:5)
    real(dp),intent(in):: x,mu
    integer,intent(in)::hadron

    if(mu<muTMM_min) then
        write(*,*) ErrorString("ERROR in KT-moment computation. mu<mu_min:",moduleName),muTMM_min
        error stop
    end if

    SiversTMDPDF_TMM_X=HankelbyOGATA%Moment_X(F,mu)

    !!! forcefully set =0 below threshold
    if(mu<mBOTTOM) then
    SiversTMDPDF_TMM_X(5)=0._dp
    SiversTMDPDF_TMM_X(-5)=0._dp
    end if
    if(mu<mCHARM) then
    SiversTMDPDF_TMM_X(4)=0._dp
    SiversTMDPDF_TMM_X(-4)=0._dp
    end if
contains
    function F(b)
    real(dp),dimension(-5:5)::F
    real(dp),intent(in)::b
    F=TMD_opt(x,b,hadron)
    end function F

end function SiversTMDPDF_TMM_X
end module SiversTMDPDF
