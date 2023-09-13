!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.00
!
!    Evaluation of the convolution integral for DY TMD cross-section with KPC
!
!    if you use this module please, quote 2307.13054
!
!    ver 3.0: created (AV, 07.09.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDF_KPC_DY
use aTMDe_Numerics
use IntegrationRoutines
use IO_functions
use EWinput
use TMDs_inKT

implicit none

private

character (len=12),parameter :: moduleName="TMDin-KPC-DY"
character (len=5),parameter :: version="v3.00"
!Last appropriate verion of constants-file
integer,parameter::inputver=30

integer::outputLevel=2
!! variable that count number of WRNING mesagges. In order not to spam too much
integer::messageTrigger=6
integer::messageCounter=0
logical::started=.false.
!! flag for loss of convergence
logical:: convergenceLost=.false.

!! tolerances for integration and general
real(dp)::toleranceGEN
real(dp)::toleranceINT

!increment counters
integer::GlobalCounter=0 !!!total counter of calls of TMD pairs
integer::LocalCounter=0 !!!counter of calls of TMD pairs within the current integrand

! parameter of TMD proportionality
real(dp)::M2=1._dp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public declarations
public::TMDF_KPC_DY_IsInitialized,TMDF_KPC_DY_Initialize
public::KPC_DYconv

contains

INCLUDE 'Code/TMDF_KPC/DY_TMDpairs.f90'
INCLUDE 'Code/TMDF_KPC/DY_KERNELpairs.f90'

function TMDF_KPC_DY_IsInitialized()
    logical::TMDF_KPC_DY_IsInitialized
    TMDF_KPC_DY_IsInitialized=started
end function TMDF_KPC_DY_IsInitialized

   !! Initialization of the package
subroutine TMDF_KPC_DY_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path
    logical::initRequired
    integer::FILEver
    real(dp)::dummy

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
        write(*,*) '             Update the const-file with artemide.setup'
        write(*,*) '  '
        CLOSE (51, STATUS='KEEP')
        stop
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger

    !!! mass parameter
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p2  ')
    read(51,*) dummy
    M2=dummy**2

    call MoveTO(51,'*14   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        CLOSE (51, STATUS='KEEP')
        return
    end if
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) toleranceGEN
    call MoveTO(51,'*p2  ')
    read(51,*) toleranceINT

    if(outputLevel>2) write(*,'(A,ES8.2)') ' | tolerance general    : ',toleranceGEN
    if(outputLevel>2) write(*,'(A,ES8.2)') ' | tolerance integral   : ',toleranceINT

    CLOSE (51, STATUS='KEEP')

    convergenceLost=.false.
    GlobalCounter=0
    LocalCounter=0
    messageCounter=0

    if(.not.TMDs_inKT_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing TMDs (from ',moduleName,')'
        if(present(prefix)) then
            call TMDs_inKT_Initialize(file,prefix)
        else
            call TMDs_inKT_Initialize(file)
        end if
    end if

    if(.not.EWinput_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing EWinput (from ',moduleName,')'
        if(present(prefix)) then
            call EWinput_Initialize(file,prefix)
        else
            call EWinput_Initialize(file)
        end if
    end if


    started=.true.
    if(outputLevel>0) write(*,*) color('----- arTeMiDe.TMDF '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '

end subroutine TMDF_KPC_DY_Initialize


!!!--------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Function that computes the integral for KPC convolution in DY
!!! proc1 = (int,int,int) is the process def for TMD*TMD
!!! proc2 = int is the process definition for the integral kernel
!!! Q, qT, x1,x2, mu are usual DY variables
!!! THIS IS A SYMMETRIC VERSION (i.e. it should contain only cos(theta)
function KPC_DYconv(Q,qT,x1,x2,mu,proc1,proc2)
    real(dp),intent(in)::Q,qT,x1,x2,mu
    integer,intent(in),dimension(1:3)::proc1
    integer,intent(in)::proc2
    real(dp)::KPC_DYconv

    real(dp)::tau2,deltaT

    LocalCounter=0

    tau2=Q**2+qT**2
    deltaT=qT**2/tau2

    KPC_DYconv=Integrate_GK(Integrand_forTheta,0._dp,pi,toleranceINT)
    write(*,*) "LC=",LocalCounter

contains

    function Integrand_forTheta(theta)
    real(dp)::Integrand_forTheta
    real(dp),intent(in)::theta
    real(dp)::cT

    cT=cos(theta)

    Integrand_forTheta=INT_overALPHA(Q,tau2,deltaT,x1,x2,mu,proc1,proc2,cT)

end function Integrand_forTheta

end function KPC_DYconv

function INT_overALPHA(Q,tau2,deltaT,x1,x2,mu,proc1,proc2,cT)
    real(dp),intent(in)::Q,tau2,deltaT,x1,x2,mu,cT
    integer,intent(in),dimension(1:3)::proc1
    integer,intent(in)::proc2
    real(dp)::INT_overALPHA

    INT_overALPHA=Integrate_GK(Integrand_forALPHA,0._dp,piHalf,toleranceINT)

contains

function Integrand_forALPHA(alpha)
    real(dp)::Integrand_forALPHA
    real(dp),intent(in)::alpha
    real(dp)::S,Lam,xi1,xi2,K1,K2,sinA,Q2

    sinA=cos(alpha)
    S=Sqrt(deltaT*sinA)*cT
    Lam=(1-deltaT)*(1-sinA)

    xi1=x1/2*(1+S+sqrt(Lam))
    xi2=x2/2*(1-S+sqrt(Lam))
    K1=tau2/4*((1+S)**2-Lam)
    K2=tau2/4*((1-S)**2-Lam)

    Q2=Q**2

    Integrand_forALPHA=DY_TMD_pair(Q,xi1,xi2,k1,k2,mu,proc1)*DY_KERNEL_pair(Q2,tau2-Q2,x1,x2,xi1,xi2,k1,k2,cT,sinA,proc2)/2


end function Integrand_forALPHA

end function INT_overALPHA

!!!--------------------------------------------------------------------------------------------------

end module TMDF_KPC_DY
