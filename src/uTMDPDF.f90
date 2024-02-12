!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	Evaluation of the unpolarized TMD PDF at low normalization point in zeta-prescription.
!	
!	if you use this module please, quote 1706.01473
!
!	18.08.2023  Implementation in ver.3.0
!
!				A.Vladimirov (18.08.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module uTMDPDF
use aTMDe_Numerics
use IO_functions
use QCDinput
use TMDR
use uTMDPDF_OPE
use uTMDPDF_model

implicit none
!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v3.00"
character (len=7),parameter :: moduleName="uTMDPDF"
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

!!!------------------------------ General parameters----------------------------------------------
logical::includeGluon=.false.   !! gluons included/non-included
integer::numOfHadrons=1         !! total number of hadrons to compute
real(dp)::TMDmass=1._dp         !! mass parameter used as mass-scale

!!!------------------------------ Parameters of transform to KT-space -------------------------------------------

integer,parameter::TMDtypeN=0 !!!!! this is the order of Bessel-transform (IT IS STRICT FOR TMD)
real(dp)::kT_FREEZE=0.0001_dp  !!!!! parameter of freezing the low-kT-value

!----Ogata Tables---
integer,parameter::Nmax=1000
INCLUDE 'Tables/BesselZero1000.f90'

!!!!! I split the qT over runs qT<qTSegmentationBoundary
!!!!! In each segment I have the ogata quadrature with h=hOGATA*hSegmentationWeight
!!!!! It helps to convergen integrals, since h(optimal) ~ qT
integer,parameter::hSegmentationNumber=7
real(dp),dimension(1:hSegmentationNumber),parameter::hSegmentationWeight=(/0.0001d0,0.001d0,0.01d0,1d0,2d0,5d0,10d0/)
real(dp),dimension(1:hSegmentationNumber),parameter::qTSegmentationBoundary=(/0.001d0,0.01d0,0.1d0,10d0,50d0,100d0,200d0/)

real(dp)::hOGATA,toleranceOGATA
!!!weights of ogata quadrature
real(dp),dimension(1:hSegmentationNumber,1:Nmax)::ww
!!!nodes of ogata quadrature
real(dp),dimension(1:hSegmentationNumber,1:Nmax)::bb

!!!------------------------------ Parameters of transform to TMM -------------------------------------------

real(dp)::muTMM_min=0.8_dp  !!!!! minimal mu

!!!!! I split the qT over runs qT<qTSegmentationBoundary
!!!!! For TMM this split is the same as for inKT

real(dp)::hOGATA_TMM,toleranceOGATA_TMM
!!!weights of ogata quadrature
real(dp),dimension(1:hSegmentationNumber,0:3,1:Nmax)::ww_TMM
!!!nodes of ogata quadrature
real(dp),dimension(1:hSegmentationNumber,0:3,1:Nmax)::bb_TMM

!!-----------------------------------------------Public interface---------------------------------------------------

public::uTMDPDF_Initialize,uTMDPDF_IsInitialized,uTMDPDF_SetScaleVariation,uTMDPDF_SetPDFreplica
public::uTMDPDF_SetLambdaNP,uTMDPDF_CurrentLambdaNP
public::uTMDPDF_lowScale5
public::uTMDPDF_inB,uTMDPDF_inKT,uTMDPDF_TMM_G,uTMDPDF_TMM_X,uPDF_uPDF

interface uTMDPDF_inB
    module procedure TMD_opt,TMD_ev
end interface

interface uTMDPDF_inKT
    module procedure Fourier_opt,Fourier_ev
end interface

contains

INCLUDE 'Code/KTspace/Fourier.f90'
INCLUDE 'Code/KTspace/Moment.f90'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Interface subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function uTMDPDF_IsInitialized()
    logical::uTMDPDF_IsInitialized
    uTMDPDF_IsInitialized=started
end function uTMDPDF_IsInitialized

!! Initialization of the package
subroutine uTMDPDF_Initialize(file,prefix)
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
        stop
    end if

    call MoveTO(51,'*4   ')
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


    if(outputLevel>2 .and. includeGluon) write(*,'(A)') ' ... gluons are included'
    if(outputLevel>2 .and. .not.includeGluon) write(*,'(A)') ' ... gluons are not included'
    if(outputLevel>2) write(*,'(A,I3)') ' Number of hadrons to be considered =',numOfHadrons
    if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',lambdaNPlength
    if(outputLevel>2) write(*,'(A,F12.2)') ' Absolute maximum b      =',BMAX_ABS

    allocate(lambdaNP(1:lambdaNPlength))

    call PrepareTables()

    if(.not.TMDR_IsInitialized()) then
        if(outputLevel>2) write(*,*) '.. initializing TMDR (from ',moduleName,')'
        if(present(prefix)) then
            call TMDR_Initialize(file,prefix)
        else
            call TMDR_Initialize(file)
        end if
    end if

    if(.not.uTMDPDF_OPE_IsInitialized()) then
        if(outputLevel>2) write(*,*) '.. initializing uTMDPDF_OPE (from ',moduleName,')'
        if(present(prefix)) then
            call uTMDPDF_OPE_Initialize(file,prefix)
        else
            call uTMDPDF_OPE_Initialize(file)
        end if
    end if

    call ModelInitialization(lambdaNPlength)
    if(outputLevel>0) write(*,*) color('----- arTeMiDe.uTMDPDF_model : .... initialized',c_green)

    started=.true.
    messageCounter=0

    if(outputLevel>0) write(*,*) color('----- arTeMiDe.uTMDPDF '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '
end subroutine uTMDPDF_Initialize

!!!!!!!!!! ------------------------ SUPPORINTG ROUTINES --------------------------------------
!!! update PDF replica
subroutine uTMDPDF_SetPDFreplica(rep,hadron)
    integer,intent(in):: rep,hadron

    call uTMDPDF_OPE_SetPDFreplica(rep,hadron)
end subroutine uTMDPDF_SetPDFreplica

!!!! this routine set the variations of scales
subroutine uTMDPDF_SetScaleVariation(c4_in)
    real(dp),intent(in)::c4_in
    call uTMDPDF_OPE_SetScaleVariation(c4_in)
end subroutine uTMDPDF_SetScaleVariation

!!!Sets the non-pertrubative parameters lambda
!!! carries additionl option to build the grid
!!! if need to build grid, specify the gluon required directive.
subroutine uTMDPDF_SetLambdaNP(lambdaIN)
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
end subroutine uTMDPDF_SetLambdaNP

!!! returns current value of NP parameters
function uTMDPDF_CurrentLambdaNP()
    real(dp),dimension(1:lambdaNPlength)::uTMDPDF_CurrentLambdaNP
    uTMDPDF_CurrentLambdaNP=lambdaNP
end function uTMDPDF_CurrentLambdaNP

!!!!!!!--------------------------- DEFINING ROUTINES ------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!======TO REMOVE
!!!!!!! the function that actually returns the uTMDPDF!
function uTMDPDF_lowScale5(x,bT,hadron)
  real(dp),dimension(-5:5)::uTMDPDF_lowScale5
  real(dp),intent(in) :: x, bT
  integer,intent(in)::hadron

  !!! test boundaries
    if(x>1d0) then
        call Warning_Raise('Called x>1 (return 0). x='//numToStr(x),messageCounter,messageTrigger,moduleName)
        uTMDPDF_lowScale5=0._dp
        return
     else if(x==1.d0) then !!! funny but sometimes FORTRAN can compare real numbers exactly
        uTMDPDF_lowScale5=0._dp
        return
    else if(bT>BMAX_ABS) then
        uTMDPDF_lowScale5=0._dp
        return
    else if(x<1d-12) then
        write(*,*) ErrorString('Called x<0. x='//numToStr(x)//' . Evaluation STOP',moduleName)
        stop
    else if(bT<0d0) then
        write(*,*) ErrorString('Called b<0. b='//numToStr(bT)//' . Evaluation STOP',moduleName)
        stop
    end if

    uTMDPDF_lowScale5=uTMDPDF_OPE_convolution(x,bT,abs(hadron))*FNP(x,bT,abs(hadron),lambdaNP)

    if(hadron<0) uTMDPDF_lowScale5=uTMDPDF_lowScale5(5:-5:-1)

end function uTMDPDF_lowScale5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!======TO REMOVE

!!!!! the names are neutral because these procedures are feed to Fourier transform. And others universal sub programs.

!!!!!!! the function that actually returns the uTMDPDF optimal value
function TMD_opt(x,bT,hadron)
  real(dp),dimension(-5:5)::TMD_opt
  real(dp),intent(in) :: x, bT
  integer,intent(in)::hadron

  !!! test boundaries
    if(x>1d0) then
        call Warning_Raise('Called x>1 (return 0). x='//numToStr(x),messageCounter,messageTrigger,moduleName)
        TMD_opt=0._dp
        return
     else if(x==1.d0) then !!! funny but sometimes FORTRAN can compare real numbers exactly
        TMD_opt=0._dp
        return
    else if(bT>BMAX_ABS) then
        TMD_opt=0._dp
        return
    else if(x<1d-12) then
        write(*,*) ErrorString('Called x<0. x='//numToStr(x)//' . Evaluation STOP',moduleName)
        stop
    else if(bT<0d0) then
        write(*,*) ErrorString('Called b<0. b='//numToStr(bT)//' . Evaluation STOP',moduleName)
        stop
    end if

    TMD_opt=uTMDPDF_OPE_convolution(x,bT,abs(hadron))*FNP(x,bT,abs(hadron),lambdaNP)

    if(hadron<0) TMD_opt=TMD_opt(5:-5:-1)

end function TMD_opt
!
!!!!!!!! the function that actually returns the uTMDPDF evolved to (mu,zeta) value
function TMD_ev(x,bt,muf,zetaf,hadron)
    real(dp)::TMD_ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: Rkernel,RkernelG

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
    TMD_ev(5)=0_dp
    TMD_ev(-5)=0_dp
    end if
    if(muf<mCHARM) then
    TMD_ev(4)=0_dp
    TMD_ev(-4)=0_dp
    end if

end function TMD_ev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PRODUCTS FOR DY!!!!!!!!!!!!!!!!!!!!!!
!!! Product of quark*antiquark uTMDPDFs. Slightly faster then just product
! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
function uPDF_uPDF(x1,x2,bt,muf,zetaf,hadron1,hadron2)
    real(dp),dimension(-5:5)::uPDF_uPDF
    real(dp),intent(in):: x1,x2,bt,muf,zetaf
    integer,intent(in)::hadron1,hadron2
    real(dp):: Rkernel, RkernelG
    real(dp),dimension(-5:5)::tmd1,tmd2


    tmd1=uTMDPDF_lowScale5(x1,bT,hadron1)
    tmd2=uTMDPDF_lowScale5(x2,bT,hadron2)

    !!!! both hadrons, so it is q*barq product
    if(sign(1,hadron1)==sign(1,hadron1)) then
        uPDF_uPDF=tmd1*(tmd2(5:-5:-1))
    else !!!!! one hadron is anti-hadron, so it is q*q product
        uPDF_uPDF=tmd1*tmd2
    end if

    if(includeGluon) then
        !!! the evolution factor is squared because There is a product of two TMDs
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)**2
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)**2

        uPDF_uPDF=uPDF_uPDF*&
            (/Rkernel,Rkernel,Rkernel,Rkernel,Rkernel,RkernelG,Rkernel,Rkernel,Rkernel,Rkernel,Rkernel/)

    else
        !!! the evolution factor is squared because There is a product of two TMDs
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)**2
        uPDF_uPDF=Rkernel*uPDF_uPDF
    end if

    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
        uPDF_uPDF(5)=0_dp
        uPDF_uPDF(-5)=0_dp
    end if
    if(muf<mCHARM) then
        uPDF_uPDF(4)=0_dp
        uPDF_uPDF(-4)=0_dp
    end if
end function uPDF_uPDF

!!!!!!!! TMM G_{n,n} at (x,mu)
function uTMDPDF_TMM_G(x,mu,hadron)
    real(dp)::uTMDPDF_TMM_G(-5:5)
    real(dp),intent(in):: x,mu
    integer,intent(in)::hadron

    uTMDPDF_TMM_G=Moment_G(x,mu,hadron)

    !!! forcefully set =0 below threshold
    if(mu<mBOTTOM) then
    uTMDPDF_TMM_G(5)=0_dp
    uTMDPDF_TMM_G(-5)=0_dp
    end if
    if(mu<mCHARM) then
    uTMDPDF_TMM_G(4)=0_dp
    uTMDPDF_TMM_G(-4)=0_dp
    end if

end function uTMDPDF_TMM_G

!!!!!!!! TMM G_{n+1,n} at (x,mu)
function uTMDPDF_TMM_X(x,mu,hadron)
    real(dp)::uTMDPDF_TMM_X(-5:5)
    real(dp),intent(in):: x,mu
    integer,intent(in)::hadron

    uTMDPDF_TMM_X=Moment_X(x,mu,hadron)

    !!! forcefully set =0 below threshold
    if(mu<mBOTTOM) then
    uTMDPDF_TMM_X(5)=0_dp
    uTMDPDF_TMM_X(-5)=0_dp
    end if
    if(mu<mCHARM) then
    uTMDPDF_TMM_X(4)=0_dp
    uTMDPDF_TMM_X(-4)=0_dp
    end if

end function uTMDPDF_TMM_X

end module uTMDPDF
