!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 1.4
!
!    Evaluation of the TMDs in KT space
!
!    if you use this module please, quote 1902.08474
!
!    ver 1.0: release (AV, 23.12.2018)
!    ver 2.00: release (AV, 29.03.2019)
!
!                A.Vladimirov (23.12.2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDs_inKT
use aTMDe_Numerics
use IO_functions
use TMDs
implicit none

private

character (len=10),parameter :: moduleName="TMDs-inKT"
character (len=5),parameter :: version="v3.00"
!Last appropriate verion of constants-file
integer,parameter::inputver=1

!------------------------------------------Tables-----------------------------------------------------------------------
integer,parameter::Nmax=1000
INCLUDE 'Tables/BesselZero1000.f90'

logical:: convergenceLost=.false.

!!!!! I split the qT over runs qT<qTSegmentationBoundary
!!!!! In each segment I have the ogata quadrature with h=hOGATA*hSegmentationWeight
!!!!! It helps to convergen integrals, since h(optimal) ~ qT
integer,parameter::hSegmentationNumber=7
real(dp),dimension(1:hSegmentationNumber),parameter::hSegmentationWeight=(/0.0001d0,0.001d0,0.01d0,1d0,2d0,5d0,10d0/)
real(dp),dimension(1:hSegmentationNumber),parameter::qTSegmentationBoundary=(/0.001d0,0.01d0,0.1d0,10d0,50d0,100d0,200d0/)

real(dp)::hOGATA,tolerance
!!!weights of ogata quadrature
real(dp),dimension(1:hSegmentationNumber,0:3,1:Nmax)::ww
!!!nodes of ogata quadrature
real(dp),dimension(1:hSegmentationNumber,0:3,1:Nmax)::bb

integer::GlobalCounter
integer::CallCounter
integer::MaxCounter
integer::messageCounter
 
!!!! general interface to the optimal TMD PDF
abstract interface
    function tmd_opt(x_f,b_f,h_f)
        use aTMDe_Numerics, only: dp
        real(dp),intent(in)::x_f,b_f
        integer,intent(in)::h_f
        real(dp),dimension(-5:5)::tmd_opt
    end function tmd_opt
end interface

!!!! general interface to the TMD PDF
abstract interface
    function tmd_ev(x_f,b_f,mu,zeta,h_f)
        use aTMDe_Numerics, only: dp
        real(dp),intent(in)::x_f,b_f,mu,zeta
        integer,intent(in)::h_f
        real(dp),dimension(-5:5)::tmd_ev
    end function tmd_ev
end interface
!------------------------------------------Physical and mathematical constants------------------------------------------

!------------------------------------------Grid parameters -------------------------------------------------------------

!!!! size in x
integer::Nx=50
!!!! minimum x
real(dp)::xMin=0.00001_dp
!!!! parameters of grid
real(dp)::DeltaX
real(dp)::parX=2._dp

!!!! size in Q
integer::NQ=25
!!!! minimum & maximum Q
real(dp)::QMin=1._dp
real(dp)::QMax=200._dp
!!!! parameters of grid
real(dp)::DeltaQ

!!!! size in kT (at Qmax)
integer::NK=50
!!!! minimum Q
real(dp)::KMin=0.0001_dp
!!!! Parameter of defining maximum value of K (=park*Q)
real(dp)::parK=2.5_dp
!!!! parameters of grid
real(dp)::DeltaK

!!!! grid matrices
real(dp),allocatable::grid_uTMDPDF(:,:,:,:,:)
real(dp),allocatable::grid_uTMDFF(:,:,:,:,:)
real(dp),allocatable::grid_lpTMDPDF(:,:,:,:,:)
real(dp),allocatable::grid_SiversTMDPDF(:,:,:,:,:)
real(dp),allocatable::grid_wgtTMDPDF(:,:,:,:,:)
real(dp),allocatable::grid_BoerMuldersTMDPDF(:,:,:,:,:)

!!!! logical statements of including grids
logical::include_grid_uTMDPDF,include_grid_uTMDFF,include_grid_lpTMDPDF,include_grid_SiversTMDPDF,&
        include_grid_wgtTMDPDF,include_grid_BoerMuldersTMDPDF
!!!! number of hadron in the grid for each type of TMD
integer::numH_uTMDPDF,numH_uTMDFF,numH_lpTMDPDF,numH_SiversTMDPDF,&
        numH_wgtTMDPDF,numH_BoerMuldersTMDPDF
!!!! logical statements of including gluons into the grid
logical::include_gluon_uTMDPDF,include_gluon_uTMDFF,include_gluon_lpTMDPDF,include_gluon_SiversTMDPDF,&
        include_gluon_wgtTMDPDF,include_gluon_BoerMuldersTMDPDF

!!!! logical statements that grids are ready
logical::IsGridReady_uTMDPDF,IsGridReady_uTMDFF,IsGridReady_lpTMDPDF,IsGridReady_SiversTMDPDF,&
        IsGridReady_wgtTMDPDF,IsGridReady_BoerMuldersTMDPDF

!------------------------------------------Working variables------------------------------------------------------------
  
logical::started=.false.

integer::outputLevel=2
integer::messageTrigger=5

public::TMDs_inKT_Initialize,TMDs_inKT_ShowStatistic,TMDs_inKT_IsInitialized,TMDs_inKT_ResetCounters
public::Moment_G,Moment_X
    
real(dp),dimension(-5:5),public::uTMDPDF_kT_50,uTMDPDF_kT_5
real(dp),dimension(-5:5),public::uTMDFF_kT_5,uTMDFF_kT_50
real(dp),dimension(-5:5),public::lpTMDPDF_kT_50
real(dp),dimension(-5:5),public::SiversTMDPDF_kT_5,SiversTMDPDF_kT_50
real(dp),dimension(-5:5),public::wgtTMDPDF_kT_5,wgtTMDPDF_kT_50
real(dp),dimension(-5:5),public::BoerMuldersTMDPDF_kT_5,BoerMuldersTMDPDF_kT_50
public::testTMD_kT

interface uTMDPDF_kT_5
    module procedure uTMDPDF_kT_5_Ev,uTMDPDF_kT_5_optimal
end interface

interface uTMDPDF_kT_50
    module procedure uTMDPDF_kT_50_Ev,uTMDPDF_kT_50_optimal
end interface

interface uTMDFF_kT_5
    module procedure uTMDFF_kT_5_Ev,uTMDFF_kT_5_optimal
end interface

interface uTMDFF_kT_50
    module procedure uTMDFF_kT_50_Ev,uTMDFF_kT_50_optimal
end interface

interface lpTMDPDF_kT_50
    module procedure lpTMDPDF_kT_50_Ev,lpTMDPDF_kT_50_optimal
end interface

interface SiversTMDPDF_kT_5
    module procedure SiversTMDPDF_kT_5_Ev,SiversTMDPDF_kT_5_optimal
end interface

interface SiversTMDPDF_kT_50
    module procedure SiversTMDPDF_kT_50_Ev,SiversTMDPDF_kT_50_optimal
end interface

interface wgtTMDPDF_kT_5
    module procedure wgtTMDPDF_kT_5_Ev,wgtTMDPDF_kT_5_optimal
end interface

interface wgtTMDPDF_kT_50
    module procedure wgtTMDPDF_kT_50_Ev,wgtTMDPDF_kT_50_optimal
end interface

interface BoerMuldersTMDPDF_kT_5
    module procedure BoerMuldersTMDPDF_kT_5_Ev,BoerMuldersTMDPDF_kT_5_optimal
end interface

interface BoerMuldersTMDPDF_kT_50
    module procedure BoerMuldersTMDPDF_kT_50_Ev,BoerMuldersTMDPDF_kT_50_optimal
end interface

contains 

INCLUDE 'Code/TMDs_inKT/grid_inKT.f90'

function TMDs_inKT_IsInitialized()
    logical::TMDs_inKT_IsInitialized
    TMDs_inKT_IsInitialized=started
end function TMDs_inKT_IsInitialized

   !! Initialization of the package
subroutine TMDs_inKT_Initialize(file,prefix)
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
        write(*,*) '             Update the const-file with artemide.setup'
        write(*,*) '  '
        stop
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger

    call MoveTO(51,'*8   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        return
    end if

    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) hOGATA
        
    if(outputLevel>2) write(*,'(A,ES8.2)') ' | h for Ogata quadrature    : ',hOGATA
    if(outputLevel>2) write(*,'(A,ES8.2)') ' | tolerance            : ',tolerance

    !!!! grid parameters
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) xMin
    call MoveTO(51,'*p2  ')
    read(51,*) parX
    call MoveTO(51,'*p3  ')
    read(51,*) KMin
    call MoveTO(51,'*p4  ')
    read(51,*) parK
    call MoveTO(51,'*p5  ')
    read(51,*) QMin
    call MoveTO(51,'*p6  ')
    read(51,*) QMax
    call MoveTO(51,'*p7  ')
    read(51,*) Nx
    call MoveTO(51,'*p8  ')
    read(51,*) NK
    call MoveTO(51,'*p9  ')
    read(51,*) NQ

    !!!! Inclusion of grid for each type of TMD
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_grid_uTMDPDF
    call MoveTO(51,'*p2  ')
    read(51,*) numH_uTMDPDF
    call MoveTO(51,'*p3  ')
    read(51,*) include_gluon_uTMDPDF
    call MoveTO(51,'*p4  ')
    read(51,*) include_grid_uTMDFF
    call MoveTO(51,'*p5  ')
    read(51,*) numH_uTMDFF
    call MoveTO(51,'*p6  ')
    read(51,*) include_gluon_uTMDFF
    call MoveTO(51,'*p7  ')
    read(51,*) include_grid_lpTMDPDF
    call MoveTO(51,'*p8  ')
    read(51,*) numH_lpTMDPDF
    call MoveTO(51,'*p9  ')
    read(51,*) include_gluon_lpTMDPDF
    call MoveTO(51,'*p10 ')
    read(51,*) include_grid_SiversTMDPDF
    call MoveTO(51,'*p11 ')
    read(51,*) numH_SiversTMDPDF
    call MoveTO(51,'*p12 ')
    read(51,*) include_gluon_SiversTMDPDF
    call MoveTO(51,'*p13 ')
    read(51,*) include_grid_wgtTMDPDF
    call MoveTO(51,'*p14 ')
    read(51,*) numH_wgtTMDPDF
    call MoveTO(51,'*p15 ')
    read(51,*) include_gluon_wgtTMDPDF
    call MoveTO(51,'*p16 ')
    read(51,*) include_grid_BoerMuldersTMDPDF
    call MoveTO(51,'*p17 ')
    read(51,*) numH_BoerMuldersTMDPDF
    call MoveTO(51,'*p18 ')
    read(51,*) include_gluon_BoerMuldersTMDPDF
        
    CLOSE (51, STATUS='KEEP') 
        
    if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs-inKT: preparing Ogata tables'
    call PrepareTables()
    if(outputLevel>2) write(*,'(A,I4)') ' | Maximum number of nodes    :',Nmax
    if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs-inKT: Ogata tables prepared'
    
    GlobalCounter=0
    CallCounter=0
    MaxCounter=0
    convergenceLost=.false.
        
    if(.not.TMDs_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing TMDs (from ',moduleName,')'
        if(present(prefix)) then
            call TMDs_Initialize(file,prefix)
        else
            call TMDs_Initialize(file)
        end if
    end if

    IsGridReady_uTMDPDF=.false.
    IsGridReady_uTMDFF=.false.
    IsGridReady_lpTMDPDF=.false.
    IsGridReady_SiversTMDPDF=.false.
    IsGridReady_wgtTMDPDF=.false.
    IsGridReady_BoerMuldersTMDPDF=.false.

    !!!!!!!!!!!!!!!!! grid preparation
    if(include_grid_uTMDPDF) then
        allocate(grid_uTMDPDF(0:NX,0:NK,0:NQ,-5:5,1:numH_uTMDPDF))
        if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs-inKT: grid for uTMDPDF to be prepared'
    else
        if(outputLevel>2) write(*,*) 'arTeMiDe.TMDs-inKT: grid for uTMDPDF to be NOT prepared'
    end if

    if(include_grid_uTMDFF) then
        allocate(grid_uTMDFF(0:NX,0:NK,0:NQ,-5:5,1:numH_uTMDFF))
        if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs-inKT: grid for uTMDFF to be prepared'
    else
        if(outputLevel>2) write(*,*) 'arTeMiDe.TMDs-inKT: grid for uTMDFF to be NOT prepared'
    end if

    if(include_grid_lpTMDPDF) then
        allocate(grid_lpTMDPDF(0:NX,0:NK,0:NQ,-5:5,1:numH_lpTMDPDF))
        if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs-inKT: grid for lpTMDPDF to be prepared'
    else
        if(outputLevel>2) write(*,*) 'arTeMiDe.TMDs-inKT: grid for lpTMDPDF to be NOT prepared'
    end if

    if(include_grid_SiversTMDPDF) then
        allocate(grid_SiversTMDPDF(0:NX,0:NK,0:NQ,-5:5,1:numH_SiversTMDPDF))
        if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs-inKT: grid for SiversTMDPDF to be prepared'
    else
        if(outputLevel>2) write(*,*) 'arTeMiDe.TMDs-inKT: grid for SiversTMDPDF to be NOT prepared'
    end if

    if(include_grid_wgtTMDPDF) then
        allocate(grid_wgtTMDPDF(0:NX,0:NK,0:NQ,-5:5,1:numH_wgtTMDPDF))
        if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs-inKT: grid for wgtTMDPDF to be prepared'
    else
        if(outputLevel>2) write(*,*) 'arTeMiDe.TMDs-inKT: grid for wgtTMDPDF to be NOT prepared'
    end if

    if(include_grid_BoerMuldersTMDPDF) then
        allocate(grid_BoerMuldersTMDPDF(0:NX,0:NK,0:NQ,-5:5,1:numH_BoerMuldersTMDPDF))
        if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs-inKT: grid for BoerMuldersTMDPDF to be prepared'
    else
        if(outputLevel>2) write(*,*) 'arTeMiDe.TMDs-inKT: grid for BoerMuldersTMDPDF to be NOT prepared'
    end if
        
    started=.true.
    if(outputLevel>0) write(*,*) color('----- arTeMiDe.TMDs-inKT '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '

end subroutine TMDs_inKT_Initialize


  
subroutine TMDs_inKT_ResetCounters()
    
    convergenceLost=.false.
    GlobalCounter=0
    CallCounter=0
    MaxCounter=0
    
end subroutine TMDs_inKT_ResetCounters
  
    !!!!!!!Functions which carry the trigger on convergences.... Its used in xSec, and probably in other places.
function TMDs_inKT_IsconvergenceLost()
  logical::TMDs_inKT_IsconvergenceLost
  
  TMDs_inKT_IsconvergenceLost=convergenceLost
end function TMDs_inKT_IsconvergenceLost
  
subroutine TMDs_inKT_convergenceISlost()  
  convergenceLost=.true.
  if(outputLevel>1) write(*,*) WarningString('convergence is lost. Counters reset.',moduleName)
  call TMDs_inKT_ResetCounters()
end subroutine TMDs_inKT_convergenceISlost
  
subroutine TMDs_inKT_ShowStatistic()
    if(convergenceLost) then
        write(*,*) '         TMDs-in-kT statistics: convergence has been lost.'
    else
        write(*,'(A,ES12.3)') 'TMDs-in-kT statistics         total calls of TMDs  :  ',Real(2*GlobalCounter)
        write(*,'(A,ES12.3)') '                              total calls of TMD_F :  ',Real(CallCounter)
        write(*,'(A,F12.3)')  '                                         avarage M :  ',Real(GlobalCounter)/Real(CallCounter)
        write(*,'(A,I12)')    '                                     maximum calls :  ',MaxCounter

    end if
end subroutine TMDs_inKT_ShowStatistic


!!!Prepare tables for Ogata quadrature with given h
!!! note that the factor 1/(2pi) is taken into ww
!!! the difference between definition in TMDF and here is 1/pi
subroutine PrepareTables()
    integer::i,k,j
    real(dp)::hS!=h*hSegmentationWeight
    real(dp)::xi

    do j=1,hSegmentationNumber
    do k=0,3
    do i=1,Nmax

    hS=hOGATA*hSegmentationWeight(j)    
    xi=JZero(k,i)

    !     ww(j,k,i)=BESSEL_JN(k,bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)&
    !         *(pi*xi*hS*Cosh(xi*hS)+Sinh(pi*Sinh(xi*hS)))/(1d0+Cosh(pi*Sinh(xi*hS)))

    !!! if we too far away in xI*hS, the double exponential grow rapidly.
    !!! and for >6, it generates term 10^{300} and exceed the presision

    if(xi*hS>6.d0) then
        bb(j,k,i)=xi*Tanh(piHalf*Sinh(xi*hS))
        ww(j,k,i)=BESSEL_JN(k,bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)/pi
        
    else
        bb(j,k,i)=xi*Tanh(piHalf*Sinh(xi*hS))
        ww(j,k,i)=BESSEL_JN(k,bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)&
        *(pi*xi*hS*Cosh(xi*hS)/(2d0*Cosh(piHalf*Sinh(xi*hS))**2)+Tanh(piHalf*Sinh(xi*hS)))/pi
    end if

    end do
    end do
    end do 
end subroutine PrepareTables
 
 !--------------------------------------INTERFACES TO TMD------------------------------------------------
 
function testTMD_kT(x,qT)
    real(dp)::testTMD_kT(-5:5)
    real(dp),intent(in)::x,qT
    testTMD_kT=Fourier(x,qT,10d0,10d0,0,1) 
end function testTMD_kT

!---------------------------------------------------uTMDPDF
function uTMDPDF_kT_5_Ev(x,qT,mu,zeta,hadron)
    real(dp)::uTMDPDF_kT_5_Ev(-5:5)
    real(dp),intent(in)::x,qT,mu,zeta
    integer,intent(in)::hadron
    uTMDPDF_kT_5_Ev=Fourier(x,qT,mu,zeta,1,hadron) 
end function uTMDPDF_kT_5_Ev

function uTMDPDF_kT_50_Ev(x,qT,mu,zeta,hadron)
    real(dp)::uTMDPDF_kT_50_Ev(-5:5)
    real(dp),intent(in)::x,qT,mu,zeta
    integer,intent(in)::hadron
    uTMDPDF_kT_50_Ev=Fourier(x,qT,mu,zeta,2,hadron) 
end function uTMDPDF_kT_50_Ev

function uTMDPDF_kT_5_optimal(x,qT,hadron)
    real(dp)::uTMDPDF_kT_5_optimal(-5:5)
    real(dp),intent(in)::x,qT
    integer,intent(in)::hadron
    uTMDPDF_kT_5_optimal=Fourier(x,qT,10d0,10d0,3,hadron) 
end function uTMDPDF_kT_5_optimal

function uTMDPDF_kT_50_optimal(x,qT,hadron)
    real(dp)::uTMDPDF_kT_50_optimal(-5:5)
    real(dp),intent(in)::x,qT
    integer,intent(in)::hadron
    uTMDPDF_kT_50_optimal=Fourier(x,qT,10d0,10d0,4,hadron) 
end function uTMDPDF_kT_50_optimal

!---------------------------------------------------uTMDFF
function uTMDFF_kT_5_Ev(x,qT,mu,zeta,hadron)
    real(dp)::uTMDFF_kT_5_Ev(-5:5)
    real(dp),intent(in)::x,qT,mu,zeta
    integer,intent(in)::hadron
    uTMDFF_kT_5_Ev=Fourier(x,qT,mu,zeta,5,hadron) 
end function uTMDFF_kT_5_Ev

function uTMDFF_kT_50_Ev(x,qT,mu,zeta,hadron)
    real(dp)::uTMDFF_kT_50_Ev(-5:5)
    real(dp),intent(in)::x,qT,mu,zeta
    integer,intent(in)::hadron
    uTMDFF_kT_50_Ev=Fourier(x,qT,mu,zeta,6,hadron) 
end function uTMDFF_kT_50_Ev

function uTMDFF_kT_5_optimal(x,qT,hadron)
    real(dp)::uTMDFF_kT_5_optimal(-5:5)
    real(dp),intent(in)::x,qT
    integer,intent(in)::hadron
    uTMDFF_kT_5_optimal=Fourier(x,qT,10d0,10d0,7,hadron) 
end function uTMDFF_kT_5_optimal

function uTMDFF_kT_50_optimal(x,qT,hadron)
    real(dp)::uTMDFF_kT_50_optimal(-5:5)
    real(dp),intent(in)::x,qT
    integer,intent(in)::hadron
    uTMDFF_kT_50_optimal=Fourier(x,qT,10d0,10d0,8,hadron) 
end function uTMDFF_kT_50_optimal


!---------------------------------------------------lpTMDPDF

function lpTMDPDF_kT_50_Ev(x,qT,mu,zeta,hadron)
    real(dp)::lpTMDPDF_kT_50_Ev(-5:5)
    real(dp),intent(in)::x,qT,mu,zeta
    integer,intent(in)::hadron
    lpTMDPDF_kT_50_Ev=Fourier(x,qT,mu,zeta,9,hadron) 
end function lpTMDPDF_kT_50_Ev

function lpTMDPDF_kT_50_optimal(x,qT,hadron)
    real(dp)::lpTMDPDF_kT_50_optimal(-5:5)
    real(dp),intent(in)::x,qT
    integer,intent(in)::hadron
    lpTMDPDF_kT_50_optimal=Fourier(x,qT,10d0,10d0,10,hadron) 
end function lpTMDPDF_kT_50_optimal

!---------------------------------------------------SiversTMDPDF
function SiversTMDPDF_kT_5_Ev(x,qT,mu,zeta,hadron)
    real(dp)::SiversTMDPDF_kT_5_Ev(-5:5)
    real(dp),intent(in)::x,qT,mu,zeta
    integer,intent(in)::hadron
    SiversTMDPDF_kT_5_Ev=Fourier(x,qT,mu,zeta,11,hadron) 
end function SiversTMDPDF_kT_5_Ev

function SiversTMDPDF_kT_50_Ev(x,qT,mu,zeta,hadron)
    real(dp)::SiversTMDPDF_kT_50_Ev(-5:5)
    real(dp),intent(in)::x,qT,mu,zeta
    integer,intent(in)::hadron
    SiversTMDPDF_kT_50_Ev=Fourier(x,qT,mu,zeta,12,hadron) 
end function SiversTMDPDF_kT_50_Ev

function SiversTMDPDF_kT_5_optimal(x,qT,hadron)
    real(dp)::SiversTMDPDF_kT_5_optimal(-5:5)
    real(dp),intent(in)::x,qT
    integer,intent(in)::hadron
    SiversTMDPDF_kT_5_optimal=Fourier(x,qT,10d0,10d0,13,hadron) 
end function SiversTMDPDF_kT_5_optimal

function SiversTMDPDF_kT_50_optimal(x,qT,hadron)
    real(dp)::SiversTMDPDF_kT_50_optimal(-5:5)
    real(dp),intent(in)::x,qT
    integer,intent(in)::hadron
    SiversTMDPDF_kT_50_optimal=Fourier(x,qT,10d0,10d0,14,hadron) 
end function SiversTMDPDF_kT_50_optimal

!---------------------------------------------------wgtTMDPDF
function wgtTMDPDF_kT_5_Ev(x,qT,mu,zeta,hadron)
    real(dp)::wgtTMDPDF_kT_5_Ev(-5:5)
    real(dp),intent(in)::x,qT,mu,zeta
    integer,intent(in)::hadron
    wgtTMDPDF_kT_5_Ev=Fourier(x,qT,mu,zeta,15,hadron) 
end function wgtTMDPDF_kT_5_Ev

function wgtTMDPDF_kT_50_Ev(x,qT,mu,zeta,hadron)
    real(dp)::wgtTMDPDF_kT_50_Ev(-5:5)
    real(dp),intent(in)::x,qT,mu,zeta
    integer,intent(in)::hadron
    wgtTMDPDF_kT_50_Ev=Fourier(x,qT,mu,zeta,16,hadron) 
end function wgtTMDPDF_kT_50_Ev

function wgtTMDPDF_kT_5_optimal(x,qT,hadron)
    real(dp)::wgtTMDPDF_kT_5_optimal(-5:5)
    real(dp),intent(in)::x,qT
    integer,intent(in)::hadron
    wgtTMDPDF_kT_5_optimal=Fourier(x,qT,10d0,10d0,17,hadron) 
end function wgtTMDPDF_kT_5_optimal

function wgtTMDPDF_kT_50_optimal(x,qT,hadron)
    real(dp)::wgtTMDPDF_kT_50_optimal(-5:5)
    real(dp),intent(in)::x,qT
    integer,intent(in)::hadron
    wgtTMDPDF_kT_50_optimal=Fourier(x,qT,10d0,10d0,18,hadron) 
end function wgtTMDPDF_kT_50_optimal

!---------------------------------------------------BoerMuldersTMDPDF
function BoerMuldersTMDPDF_kT_5_Ev(x,qT,mu,zeta,hadron)
    real(dp)::BoerMuldersTMDPDF_kT_5_Ev(-5:5)
    real(dp),intent(in)::x,qT,mu,zeta
    integer,intent(in)::hadron
    BoerMuldersTMDPDF_kT_5_Ev=Fourier(x,qT,mu,zeta,19,hadron)
end function BoerMuldersTMDPDF_kT_5_Ev

function BoerMuldersTMDPDF_kT_50_Ev(x,qT,mu,zeta,hadron)
    real(dp)::BoerMuldersTMDPDF_kT_50_Ev(-5:5)
    real(dp),intent(in)::x,qT,mu,zeta
    integer,intent(in)::hadron
    BoerMuldersTMDPDF_kT_50_Ev=Fourier(x,qT,mu,zeta,20,hadron)
end function BoerMuldersTMDPDF_kT_50_Ev

function BoerMuldersTMDPDF_kT_5_optimal(x,qT,hadron)
    real(dp)::BoerMuldersTMDPDF_kT_5_optimal(-5:5)
    real(dp),intent(in)::x,qT
    integer,intent(in)::hadron
    BoerMuldersTMDPDF_kT_5_optimal=Fourier(x,qT,10d0,10d0,21,hadron)
end function BoerMuldersTMDPDF_kT_5_optimal

function BoerMuldersTMDPDF_kT_50_optimal(x,qT,hadron)
    real(dp)::BoerMuldersTMDPDF_kT_50_optimal(-5:5)
    real(dp),intent(in)::x,qT
    integer,intent(in)::hadron
    BoerMuldersTMDPDF_kT_50_optimal=Fourier(x,qT,10d0,10d0,22,hadron)
end function BoerMuldersTMDPDF_kT_50_optimal

!------------------------------------------FOURIER--------------------------------
!!!This is the defining module function
!!! It evaluates the integral 
!!!  int_0^infty   b db/2pi  J_num(b qT) F1
!!! The function F1 is given via number.. num
!!! The order of the Bessel function is also selected according to num
function Fourier(x,qT_in,mu,zeta,num,hadron)
    real(dp),intent(in)::x,mu,zeta,qT_in
    integer,intent(in)::num,hadron
    real(dp)::integral(-5:5),eps(-5:5),qT
    real(dp)::v1(-5:5),v2(-5:5),v3(-5:5),v4(-5:5),delta(-5:5)
    logical:: partDone(-5:5)
    integer::k,n,j,Nsegment
    real(dp)::Fourier(-5:5)

    CallCounter=CallCounter+1
    integral=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)

    if(qT_in<0.0001d0) then  
        qT=0.0001d0  
    else
        qT=qT_in
    end if

    !!!in the case of lost convergence we return huge number (divergent xSec)
    if(TMDs_inKT_IsconvergenceLost()) then
        Fourier=integral+1d10
    else

    v1=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v2=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v3=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v4=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    partDone=(/.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false./)

    !!! define segment of qT
    do j=1,hSegmentationNumber
        if(qT<qTSegmentationBoundary(j)) exit
    end do
    if(j>hSegmentationNumber) then
        Nsegment=hSegmentationNumber
    else
        Nsegment=j
    end if

    !!!! select the order of Bessel function for transform
    SELECT CASE(num)
        CASE(9,10)
            n=2
        CASE(11,12,13,14,15,16,17,18)
            n=1            
        CASE DEFAULT
            n=0
    END SELECT

    do k=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
        eps=ww(Nsegment,n,k)*(bb(Nsegment,n,k)**(n+1))*Integrand(bb(Nsegment,n,k)/qT,x,mu,zeta,num,hadron)
            
        v4=v3
        v3=v2
        v2=v1
        v1=abs(eps)

        delta=(v1+v2+v3+v4)
        integral=integral+eps
            
        !!! here we check that residual term is smaller than already collected integral
        !!! also checking the zerothness of the integral. If already collected integral is null it is null
        !!! Here is potential bug. If the first 10 points give zero (whereas some later points do not), the integral will be zero
        !!! I check for each separate flavor
        do j=-5,5
            if((delta(j)<tolerance*ABS(integral(j)) .or. ABS(integral(j))<1d-32) .and. k>=10) partDone(j)=.true.
        end do
        if(partDone(-5).and.partDone(-4).and.partDone(-3).and.partDone(-2).and.partDone(-1)&
            .and.partDone(0).and.partDone(1).and.partDone(2).and.partDone(3).and.partDone(4).and.partDone(5)) exit
            
    end do

    if(k>=Nmax) then
        if(outputlevel>0) call Warning_Raise('OGATA quadrature diverge. TMD decaing too slow?',&
            messageCounter,messageTrigger,moduleName)
            if(outputlevel>2) then
            write(*,*) 'Information over the last call ----------'
            write(*,*) partDone
            write(*,*) 'bt/qT= ',bb(Nsegment,n,Nmax)/qT, 'qT=',qT, '| segmentation zone=',Nsegment,&
                ' ogata h=',hOGATA*hSegmentationWeight(Nsegment)
            write(*,*) 'W=',Integrand(bb(Nsegment,n,Nmax)/qT,x,mu,zeta,num,hadron), 'eps/integral =', eps/integral
            write(*,*) 'v1+v2+v3+v4=',v1+v2+v3+v4, '>',tolerance*(ABS(integral(1))+ABS(integral(2)))
            write(*,*) 'x=',x,'type =',num,' it is ',CallCounter,' call.'
            write(*,*) '------------------------------------------'
            end if
        call TMDs_inKT_convergenceISlost()
    end if

    !! store the maximum number of calls
    if(k>MaxCounter) MaxCounter=k-1

    !!! result is scaled by qT [because the argument of Bessel was scaled bqT-> B]
    Fourier=integral/(qT**(n+2))
    end if 
    !write(*,*) 'Last call: ',k
end function Fourier
 
function Integrand(b,x,mu,zeta,num,hadron)
    real(dp),intent(in)::b,x,mu,zeta
    integer,intent(in)::num,hadron
    real(dp)::Integrand(-5:5)

    !increment counter 
    GlobalCounter=GlobalCounter+1
    SELECT CASE(num)
        CASE(0) !!! test case
            Integrand=(/(b**3d0)*Exp(-x*b),(b**2d0)*Exp(-x*b),b*Exp(-x*b),Exp(-x*b),&
            1d0/(b**2+x**2d0)**4d0,1d0/(b**2+x**2d0)**2d0,1d0/(b**2+x**2d0),&
            (b**2d0)*Exp(-x*b*b),Exp(-x*b*b), BESSEL_J0(x/b)/b,1d0/b/)

        CASE(1) !!! uTMDPDF  quarks
            Integrand=uTMDPDF_5(x,b,mu,zeta,hadron)
            Integrand(0)=0d0

        CASE(2) !!! uTMDPDF  quarks+gluon
            Integrand=uTMDPDF_50(x,b,mu,zeta,hadron)

        CASE(3) !!! uTMDPDF  quarks OPTIMAL
            Integrand=uTMDPDF_5(x,b,hadron)
            Integrand(0)=0d0

        CASE(4) !!! uTMDPDF  quarks+gluon OPTIMAL
            Integrand=uTMDPDF_50(x,b,hadron)

        CASE(5) !!! uTMDFF  quarks
            Integrand=uTMDFF_5(x,b,mu,zeta,hadron)
            Integrand(0)=0d0

        CASE(6) !!! uTMDFF  quarks+gluon
            Integrand=uTMDFF_50(x,b,mu,zeta,hadron)

        CASE(7) !!! uTMDFF  quarks OPTIMAL
            Integrand=uTMDFF_5(x,b,hadron)
            Integrand(0)=0d0

        CASE(8) !!! uTMDFF  quarks+gluon OPTIMAL
            Integrand=uTMDFF_50(x,b,hadron)

        CASE(9) !!! lin.pol.gluon TMDPDF
            !!! minus is due to definition (see manual)
            Integrand=-lpTMDPDF_50(x,b,mu,zeta,hadron)

        CASE(10) !!! lim.pol.gluon TMDPDF OPTIMAL
            !!! minus is due to definition (see manual)
            Integrand=-lpTMDPDF_50(x,b,hadron)

        CASE(11) !!! SiversTMDPDF  quarks
            Integrand=SiversTMDPDF_5(x,b,mu,zeta,hadron)
            Integrand(0)=0d0

        CASE(12) !!! SiversuTMDPDF  quarks+gluon
            Integrand=SiversTMDPDF_50(x,b,mu,zeta,hadron)

        CASE(13) !!! SiversTMDPDF  quarks OPTIMAL
            Integrand=SiversTMDPDF_5(x,b,hadron)
            Integrand(0)=0d0

        CASE(14) !!! SiversTMDPDF  quarks+gluon OPTIMAL
            Integrand=SiversTMDPDF_50(x,b,hadron)
            
        CASE(15) !!! wgtTMDPDF  quarks
            Integrand=wgtTMDPDF_5(x,b,mu,zeta,hadron)
            Integrand(0)=0d0

        CASE(16) !!! wgtTMDPDF  quarks+gluon
            Integrand=wgtTMDPDF_50(x,b,mu,zeta,hadron)

        CASE(17) !!! wgtTMDPDF  quarks OPTIMAL
            Integrand=wgtTMDPDF_5(x,b,hadron)
            Integrand(0)=0d0

        CASE(18) !!! wgtTMDPDF  quarks+gluon OPTIMAL
            Integrand=wgtTMDPDF_50(x,b,hadron)

        CASE(19) !!! BoerMuldersTMDPDF  quarks
            Integrand=BoerMuldersTMDPDF_5(x,b,mu,zeta,hadron)
            Integrand(0)=0d0

        CASE(20) !!! BoerMuldersuTMDPDF  quarks+gluon
            Integrand=BoerMuldersTMDPDF_50(x,b,mu,zeta,hadron)

        CASE(21) !!! BoerMuldersTMDPDF  quarks OPTIMAL
            Integrand=BoerMuldersTMDPDF_5(x,b,hadron)
            Integrand(0)=0d0

        CASE(22) !!! BoerMuldersTMDPDF  quarks+gluon OPTIMAL
            Integrand=BoerMuldersTMDPDF_50(x,b,hadron)

        CASE DEFAULT
        write(*,*) ErrorString('undefined TMD: ',moduleName)
        write(*,*) color('Evaluation stop',c_red_bold)
        stop
    END SELECT

end function Integrand
 
!-----------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------functions for Computation of moments ---------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------

!------------------------------------------General Moment--------------------------------
!!! This is the general function for the computation of the Moment
!!! It evaluates the integral
!!! int_0^infty   b^n db   J_k(b mu) F[x,b,optimal]
!!! The function F is profided via directly it must be of the type uTMDPDF_lowScale5(x,bT,hadron)
function Moment_Gen(n,k,x,mu,F_opt,hadron)
    real(dp),intent(in)::x,mu
    integer,intent(in)::n,k,hadron
    procedure(tmd_opt)::F_opt
    real(dp),dimension(-5:5)::Moment_Gen(-5:5)
    real(dp),dimension(-5:5)::integral,eps,v1,v2,v3,v4,delta
    logical:: partDone(-5:5)
    integer::r,j,Nsegment

    CallCounter=CallCounter+1
    integral=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)

    if(mu<0.8d0) then
        write(*,*) ErrorString("ERROR in KT-moment computation. mu<0.8",moduleName)
        stop
    end if

    if(k>3) then
        write(*,*) ErrorString("ERROR in KT-moment computation. Called moment with J_k, k>3",moduleName)
        stop
    end if

    if(k+n<0) then
        write(*,*) ErrorString("ERROR in KT-moment computation. Integral is divergent at 0",moduleName)
        stop
    end if

    !!!in the case of lost convergence we return huge number (divergent xSec)
    if(TMDs_inKT_IsconvergenceLost()) then
        Moment_Gen=integral+1d10
    else

    v1=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v2=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v3=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v4=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    partDone=(/.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false./)

    !!! define segment of qT
    do j=1,hSegmentationNumber
        if(mu<qTSegmentationBoundary(j)) exit
    end do
    if(j>hSegmentationNumber) then
        Nsegment=hSegmentationNumber
    else
        Nsegment=j
    end if

    do r=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
        eps=ww(Nsegment,k,r)*(bb(Nsegment,k,r)**n)*F_opt(x,bb(Nsegment,k,r)/mu,hadron)

        v4=v3
        v3=v2
        v2=v1
        v1=abs(eps)

        delta=(v1+v2+v3+v4)
        integral=integral+eps

        !!! here we check that residual term is smaller than already collected integral
        !!! also checking the zerothness of the integral. If already collected integral is null it is null
        !!! Here is potential bug. If the first 10 points give zero (whereas some later points do not), the integral will be zero
        !!! I check for each separate flavor
        do j=-5,5
            if((delta(j)<tolerance*ABS(integral(j)) .or. ABS(integral(j))<1d-32) .and. r>=10) partDone(j)=.true.
        end do
        if(partDone(-5).and.partDone(-4).and.partDone(-3).and.partDone(-2).and.partDone(-1)&
            .and.partDone(0).and.partDone(1).and.partDone(2).and.partDone(3).and.partDone(4).and.partDone(5)) exit

    end do

    if(r>=Nmax) then
        if(outputlevel>0) call Warning_Raise('OGATA quadrature diverge for the Momen. TMD decaing too slow?',&
            messageCounter,messageTrigger,moduleName)
            if(outputlevel>2) then
            write(*,*) 'Information over the last call ----------'
            write(*,*) partDone
            write(*,*) 'bt/mu= ',bb(Nsegment,k,Nmax)/mu, 'qT=',mu, '| segmentation zone=',Nsegment,&
                ' ogata h=',hOGATA*hSegmentationWeight(Nsegment)
            write(*,*) 'W=',F_opt(x,bb(Nsegment,k,Nmax)/mu,hadron), 'eps/integral =', eps/integral
            write(*,*) 'v1+v2+v3+v4=',v1+v2+v3+v4, '>',tolerance*(ABS(integral(1))+ABS(integral(2)))
            write(*,*) 'x=',x,'n=',n,'type of J =',k,' it is ',CallCounter,' call.'
            write(*,*) '------------------------------------------'
            end if
        call TMDs_inKT_convergenceISlost()
    end if

    !! store the maximum number of calls
    if(r>MaxCounter) MaxCounter=r-1

    !!! result is multiplied by (2pi) [because the waights are defined for db/2pi] and scaled by mu [by definition of the moment]
    Moment_Gen=integral/(mu**(n+1))*pix2
    end if
    !write(*,*) 'Last call: ',r
end function Moment_Gen


!!!! The moment G_n, defined as
!!!! mu^{n+1}/2^n/n! int b^n J_{n+1} (b mu) F[OPTIMAL]
function Moment_G(n,x,mu,F_opt,hadron)
    real(dp),intent(in)::x,mu
    integer,intent(in)::n,hadron
    procedure(tmd_opt)::F_opt
    real(dp),dimension(-5:5)::Moment_G(-5:5)

    SELECT CASE(n)
        CASE(0)
            Moment_G=mu*Moment_Gen(0,1,x,mu,F_opt,hadron)
        CASE(1)
            Moment_G=mu**2/2*Moment_Gen(1,2,x,mu,F_opt,hadron)
        CASE(2)
            Moment_G=mu**3/8*Moment_Gen(2,3,x,mu,F_opt,hadron)
        CASE DEFAULT
            write(*,*) ErrorString("MOMENT G is defined only for n=0,1,2",moduleName)
            stop
    END SELECT

end function Moment_G

!!!! The moment G_n, defined as
!!!! mu^{n+1}/2^n/n! int b^n J_{n+1} (b mu) F[OPTIMAL]
function Moment_X(n,x,mu,F_opt,hadron)
    real(dp),intent(in)::x,mu
    integer,intent(in)::n,hadron
    procedure(tmd_opt)::F_opt
    real(dp),dimension(-5:5)::Moment_X(-5:5)

    SELECT CASE(n)
        CASE(0)
            !Moment_X=mu**3*(Moment_Gen(0,1,x,mu,F_opt,hadron)-Moment_Gen(0,3,x,mu,F_opt,hadron))/2
            Moment_X=mu**2*(mu*Moment_Gen(0,1,x,mu,F_opt,hadron)-2*Moment_Gen(-1,2,x,mu,F_opt,hadron))
        CASE(1)
            Moment_X=mu**3/2*(mu*Moment_Gen(1,2,x,mu,F_opt,hadron)-2*Moment_Gen(0,3,x,mu,F_opt,hadron))
        CASE DEFAULT
            write(*,*) ErrorString("MOMENT X is defined only for n=0,1",moduleName)
            stop
    END SELECT

end function Moment_X

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! GRIDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! make the grid
!!!!! the optional argument specifies what to do. If not present all grids (which are included) are updated
subroutine UpdateGrid(t)
character(len=3),intent(in),optional::t

if(present(t)) then

    SELECT CASE(t)
        CASE('upl')
            if(include_grid_uTMDPDF) then
                if(include_gluon_uTMDPDF) then
                    call MakeGrid(grid_uTMDPDF,2,numH_uTMDPDF)
                else
                    call MakeGrid(grid_uTMDPDF,1,numH_uTMDPDF)
                end if

                IsGridReady_uTMDPDF=.true.
            end if

        CASE('uFF')
            if(include_grid_uTMDFF) then
                if(include_gluon_uTMDFF) then
                    call MakeGrid(grid_uTMDFF,6,numH_uTMDFF)
                else
                    call MakeGrid(grid_uTMDFF,5,numH_uTMDFF)
                end if

                IsGridReady_uTMDFF=.true.
            end if

        CASE('lp ')
            if(include_grid_lpTMDPDF) then
                call MakeGrid(grid_lpTMDPDF,9,numH_lpTMDPDF)

                IsGridReady_lpTMDPDF=.true.
            end if

        CASE('Siv')
            if(include_grid_SiversTMDPDF) then
                if(include_gluon_SiversTMDPDF) then
                    call MakeGrid(grid_SiversTMDPDF,12,numH_SiversTMDPDF)
                else
                    call MakeGrid(grid_SiversTMDPDF,11,numH_SiversTMDPDF)
                end if

                IsGridReady_SiversTMDPDF=.true.
            end if

        CASE('wgt')
            if(include_grid_wgtTMDPDF) then
                if(include_gluon_wgtTMDPDF) then
                    call MakeGrid(grid_wgtTMDPDF,16,numH_wgtTMDPDF)
                else
                    call MakeGrid(grid_wgtTMDPDF,15,numH_wgtTMDPDF)
                end if

                IsGridReady_wgtTMDPDF=.true.
            end if

        CASE('BM ')
            if(include_grid_BoerMuldersTMDPDF) then
                if(include_gluon_BoerMuldersTMDPDF) then
                    call MakeGrid(grid_BoerMuldersTMDPDF,20,numH_BoerMuldersTMDPDF)
                else
                    call MakeGrid(grid_BoerMuldersTMDPDF,19,numH_BoerMuldersTMDPDF)
                end if

                IsGridReady_BoerMuldersTMDPDF=.true.
            end if

        CASE DEFAULT
            write(*,*) ErrorString("Unknown optional argument .type.",moduleName)
            stop
    END SELECT
else
    if(include_grid_uTMDPDF) then
        if(include_gluon_uTMDPDF) then
            call MakeGrid(grid_uTMDPDF,2,numH_uTMDPDF)
        else
            call MakeGrid(grid_uTMDPDF,1,numH_uTMDPDF)
        end if
    end if

    if(include_grid_uTMDFF) then
        if(include_gluon_uTMDFF) then
            call MakeGrid(grid_uTMDFF,6,numH_uTMDFF)
        else
            call MakeGrid(grid_uTMDFF,5,numH_uTMDFF)
        end if
    end if

    if(include_grid_lpTMDPDF) then
        call MakeGrid(grid_lpTMDPDF,9,numH_lpTMDPDF)
    end if

    if(include_grid_SiversTMDPDF) then
        if(include_gluon_SiversTMDPDF) then
            call MakeGrid(grid_SiversTMDPDF,12,numH_SiversTMDPDF)
        else
            call MakeGrid(grid_SiversTMDPDF,11,numH_SiversTMDPDF)
        end if
    end if

    if(include_grid_wgtTMDPDF) then
        if(include_gluon_wgtTMDPDF) then
            call MakeGrid(grid_wgtTMDPDF,16,numH_wgtTMDPDF)
        else
            call MakeGrid(grid_wgtTMDPDF,15,numH_wgtTMDPDF)
        end if
    end if

    if(include_grid_BoerMuldersTMDPDF) then
        if(include_gluon_BoerMuldersTMDPDF) then
            call MakeGrid(grid_BoerMuldersTMDPDF,20,numH_BoerMuldersTMDPDF)
        else
            call MakeGrid(grid_BoerMuldersTMDPDF,19,numH_BoerMuldersTMDPDF)
        end if
    end if
end if
end subroutine UpdateGrid

end module TMDs_inKT
