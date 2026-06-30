!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.05
!
!    This is a part of TMDX_DY complex of modules.
!    This module computes the cross-section of a Drell-Yan reaction (and related sub-processes)
!    at a single kinematical point.
!
!    if you use this module please, quote 1706.01473
!
!    ver 3.05: created from v.3.04 (AV, 29.06.2026)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDX_DY_point
use aTMDe_Numerics
use aTMDe_IO
use aTMDe_ptSpec
use TMDF
use TMDF_KPC
use LeptonCutsDY
use QCDinput
use EWinput

implicit none
private

!Current version of module
character (len=11),parameter :: moduleName="TMDX-DY-1pt"
character (len=5),parameter :: version="v3.05"
!Last appropriate version of constants-file
integer,parameter::inputver=38

integer::outputlevel=2
type(Warning_OBJ)::Warning_Handler

!!!!! numeric precisions
real(dp) :: toleranceINT=0.0001d0
real(dp) :: toleranceGEN=0.0000001d0

!!Global definitions of the factorization theorem
integer:: orderH_global
logical::useKPC
logical::usePIresum
logical:: exactX1X2    !!!.true. = utilization of exact LP value of x1 (and x2) =q^+/P+, while .false. x1=Q/sqrt(s) [only for LP factorization]
logical:: exactScales  !!!.true. = utilization of exact LP value of zeta =2q^-q^+, while .false. zeta=Q2 [only for LP factorization]

real(dp)::c2_global

integer::GlobalCounter   !!!! this counts the total number of calls of xSec in the whole aTMDe session
integer::CallCounter     !!!! this counts the total number of calls of xSec after the last reset.

real(dp)::hc2

logical::started=.false.

public:: TMDX_DY_1pt_Initialize, TMDX_DY_1pt_IsInitialized
public:: TMDX_DY_1pt_ResetCounters, TMDX_DY_1pt_SetScaleVariation
public:: xSec_DY_1pt,isConvergenceLost_inSubmodule

contains

function TMDX_DY_1pt_IsInitialized()
  logical::TMDX_DY_1pt_IsInitialized
  TMDX_DY_1pt_IsInitialized=started
end function TMDX_DY_1pt_IsInitialized

!! Initialization of the package
subroutine TMDX_DY_1pt_Initialize(file,prefix)
    character(len=*),intent(in)::file
    character(len=*),intent(in),optional::prefix
    character(len=:),allocatable::path
    logical::initRequired
    character(len=8)::orderMain
    integer::FILEver,messageTrigger

    if(started) return

    if(present(prefix)) then
    path=trim(adjustl(prefix))//trim(adjustl(file))
    else
    path=trim(adjustl(file))
    end if

    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")

    !!! Check the file version
    call MoveTO(51,'*0   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) FILEver
    if(FILEver<inputver) then
    write(*,*) 'artemide.'//trim(moduleName)//': const-file version is too old.'
    write(*,*) '             Update the const-file with artemide.setup'
    write(*,*) '  '
    CLOSE (51, STATUS='KEEP')
    error stop
    end if

    !!! Fill the message system
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger

    !!! other variables
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) hc2

    !!! go to section TMD-DY
    call MoveTO(51,'*9   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
    if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not required. '
    started=.false.
    CLOSE (51, STATUS='KEEP')
    return
    end if

    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) orderMain
    SELECT CASE(orderMain)
    CASE ("LO")
    orderH_global=0
    CASE ("NLO")
    orderH_global=1
    CASE ("NNLO")
    orderH_global=2
    CASE ("N2LO")
    orderH_global=2
    CASE ("NNNLO")
    orderH_global=3
    CASE ("N3LO") !!! same as NNNLO
    orderH_global=3
    CASE ("N4LO")
    orderH_global=4
    CASE DEFAULT
    if(outputLevel>0) write(*,*)  WarningString('try to set unknown order. Switch to NLO.',moduleName)
    orderH_global=1
    END SELECT
    if(outputLevel>1) write(*,*) '    artemide.TMDX_DY: the used order is ',trim(orderMain)

    call MoveTO(51,'*p2  ')
    read(51,*) useKPC
    if(outputLevel>1 .and. useKPC) write(*,*) color('    artemide.TMDX_DY: using TMD factorization with KPC',c_cyan)
    if(outputLevel>1 .and. useKPC) write(*,*) color('                      Please, cite [2307.13054]',c_cyan)
    if(outputLevel>1 .and. .not.(useKPC)) write(*,*) color('    artemide.TMDX_DY: using TMD factorization at LP',c_cyan)

    !! pi2 resummation
    call MoveTO(51,'*p3  ')
    read(51,*) usePIresum
    if(outputLevel>2 .and. usePIresum) write(*,*) '    artemide.TMDX_DY: pi-resummation in coef.function included.'

    !!!------ parameters of numerics
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) toleranceGEN
    call MoveTO(51,'*p2  ')
    read(51,*) toleranceINT

    if(.not.(useKPC)) then
        !!!------ parameters of LP factorization
        call MoveTO(51,'*C   ')
        !!exact values of x1x2
        call MoveTO(51,'*p1  ')
        read(51,*) exactX1X2
        if(outputLevel>2 .and. exactX1X2) &
            write(*,*) '    artemide.TMDX_DY: qT/Q corrections for x1 and x2 variables are included.'
        !!exact values for scales
        call MoveTO(51,'*p2  ')
        read(51,*) exactScales
        if(outputLevel>2 .and. exactScales) &
            write(*,*) '    artemide.TMDX_DY: qT/Q correction for scales variables are included.'

    else
        !!!------ parameters of KPC factorization
        exactX1X2=.true.
        exactScales=.false.

        call MoveTO(51,'*D   ')
        !!!!! nothing is here yet

    end if

    CLOSE (51, STATUS='KEEP')

    Warning_Handler=Warning_OBJ(moduleName=moduleName,messageCounter=0,messageTrigger=messageTrigger)

    if(.not.EWinput_IsInitialized()) then
    if(outputLevel>1) write(*,*) '.. initializing EWinput (from ',moduleName,')'
    if(present(prefix)) then
        call EWinput_Initialize(file,prefix)
    else
        call EWinput_Initialize(file)
    end if
    end if

    if(useKPC) then
    if(.not.TMDF_KPC_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing TMDF (from ',moduleName,')'
        if(present(prefix)) then
        call TMDF_KPC_Initialize(file,prefix)
        else
        call TMDF_KPC_Initialize(file)
        end if
    end if
    else
    if(.not.TMDF_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing TMDF (from ',moduleName,')'
        if(present(prefix)) then
        call TMDF_Initialize(file,prefix)
        else
        call TMDF_Initialize(file)
        end if
    end if
    end if

    !!!!! initializing Lepton Cut module
    call InitializeLeptonCutDY(toleranceINT,toleranceGEN)

    c2_global=1d0

    GlobalCounter=0
    CallCounter=0

    started=.true.

    write(*,*)  color('----- arTeMiDe.TMD_DY_1pt '//trim(version)//': .... initialized',c_green)
end subroutine TMDX_DY_1pt_Initialize

subroutine TMDX_DY_1pt_ResetCounters()
    if(outputlevel>2) then
        write(*,'(A,ES12.3)') 'TMDX DY statistics      total calls of xSec (complete session) :  ',Real(GlobalCounter)
        write(*,'(A,ES12.3)') '                        total calls of xSec (current NP-setup) :  ',Real(CallCounter)
    end if
    CallCounter=0
    call Warning_Handler%Reset()
end subroutine TMDX_DY_1pt_ResetCounters

!!!!Call this after TMD initialization but before NP, and X parameters
subroutine TMDX_DY_1pt_SetScaleVariation(c2_in)
  real(dp),intent(in)::c2_in

  if(outputLevel>1) write(*,*) 'TMDX_DY: c2 scale reset:',c2_in

  if(c2_in<0.1d0 .or. c2_in>10.d0) then
    if(outputLevel>0) write(*,*) WarningString('variation in c2 is enourmous. c2 is set to 2',moduleName)
      c2_global=2d0
  else
    c2_global=c2_in
  end if

end subroutine TMDX_DY_1pt_SetScaleVariation

!!!!! wrapper, to check the convergence in proper module.
function isConvergenceLost_inSubmodule()
logical::isConvergenceLost_inSubmodule
if(useKPC) then
  isConvergenceLost_inSubmodule=TMDF_KPC_IsconvergenceLost()
else
  isConvergenceLost_inSubmodule=TMDF_IsconvergenceLost()
end if
end function isConvergenceLost_inSubmodule


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS FOR PREFACTORS & PROCESSES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PROCESS DEFINITION:
!!! (i, h1, h2, n1) -- integers
!!! i = the integration type and global prefactor
!!! h1,h2 = the hadron types.
!!! n = the integrand combination
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! for easier read coeff-functions are split into separate file
INCLUDE '/Code/TMDX/DYcoeff-func.f90'

!!!!!----------------------------------------------------------------
!!!!!  Prefactors and cuts for LP computation
!!!!!----------------------------------------------------------------

!!!!! Prefactor 2 is (universal part) x H
!!!!! This depends only on the integration type of the process (i)
!!!!! proc1 is the first number of process array
!!!!! arguments are kinematic specification, the factorization scale, and the process
function PreFactor2(ptKinematic,scaleMu,proc1)
  type(DYpoint),intent(in)::ptKinematic
  integer,intent(in)::proc1
  real(dp),intent(in)::scaleMu
  real(dp)::PreFactor2

  real(dp)::uniPart,tauX1X2

  SELECT CASE(proc1)
  case(0)
    uniPart=1d0

  CASE(1)
    !4 pi aEm^2/3 /Nc/Q^2/s
    uniPart=pix4/9d0*(alphaEM(scaleMu)**2)/(ptKinematic%s*ptKinematic%Q2)*&
        HardCoefficientDY(scaleMu)*&
        hc2*1d9!from GeV to pb

  CASE(2)
    !! in a rare case of non-exact kinematic one should update this variable
    if(exactX1X2) then
        tauX1X2=ptKinematic%tau
    else
        tauX1X2=ptKinematic%Q/sqrt(ptKinematic%s)
    end if
    !4 pi aEm^2/3 /Nc/Q^2/s
    !! the process=2 is for the xF-integration. It has extra weight J^{-1} with J=2sqrt[(Q^2+q_T^2)/s] Cosh[y]
    !!
    uniPart=pix4/9d0*(alphaEM(scaleMu)**2)/(ptKinematic%s*ptKinematic%Q2)*&
        HardCoefficientDY(scaleMu)*&
        hc2*1d9*&!from GeV to pb
        1./(2._dp*tauX1X2*coshY(ptKinematic))  !!!! 1/J

  CASE (3) !Zboson in the narrow-width approximation
    !4 pi^2 aem/Nc/s Br(z->ee+mumu)
    uniPart=pi2x4/3d0*alphaEM(scaleMu)/ptKinematic%s*&
        HardCoefficientDY(scaleMu)*&
        hc2*1d9*&!from GeV to pb
        0.03645d0!Br from PDG, ee+mumu

  CASE (4) !Wboson in the narrow-width approximation
    !4 pi^2 aem/Nc/s Br(w->leptons)
    uniPart=pi2x4/3d0*alphaEM(scaleMu)/ptKinematic%s*&
        HardCoefficientDY(scaleMu)*&
        hc2*1d9*&!from GeV to pb
        0.1086d0!Br from PDG, ee+mumu

  CASE (5) !exclusive HIGGSboson production
    ! (2\pi) *pi Mh^2 as(mu)/36/s/vev^2 * H*cT^2
    ! (1.033)^2 is correction for mT mass in Ct at LO.
    uniPart=(1d0/18d0)*MH2*(As(c2_global*scaleMu)/VEVH)**2/ptKinematic%s*&
        HardCoefficientHIGGS(scaleMu)*(EffCouplingHFF(scaleMu)**2)*1.0677023627519822d0*&
        hc2*1d9!from GeV to pb

  CASE DEFAULT
    error stop ErrorString('unknown process p='//numToStr(proc1)//' .Evaluation stop.',moduleName)
  END SELECT

  PreFactor2=uniPart

end function PreFactor2

!!!! Fiducial cut for leptons.
!!!! this factor depends on the type of process and computed in the corresponding module
!!!! This function perform the preliminary selection depending on the process
function LeptonCutFactorLP(ptKinematic,proc1, includeCuts_in,CutParam)
  type(DYpoint),intent(in)::ptKinematic
  logical,intent(in)::includeCuts_in
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,intent(in)::proc1
  real(dp)::LeptonCutFactorLP

  !!!!! lepton-cut prefactor
  if(includeCuts_in) then
    !!! here include cuts on the lepton tensor
    LeptonCutFactorLP=CutFactor(qT_in=ptKinematic%qT,Q_in=ptKinematic%Q,y_in=ptKinematic%y,CutParameters=CutParam,Cut_Type=-2)
  else
    !!! this is uncut lepton tensor
    LeptonCutFactorLP=(1+0.5d0*(ptKinematic%qT/ptKinematic%Q)**2)
  end if

end function LeptonCutFactorLP

!!!!!----------------------------------------------------------------
!!!!!  Prefactors and cuts for KPC computation
!!!!!----------------------------------------------------------------

!!!!! Prefactor is (universal part)  x H
!!!!! This prefactor depends only on the general class of process
!!!!! proc1 is the first number of process array
!!!!! arguments are kinematic specification, the factorization scale, and the process
!!!!! Lepton cut factor is defined in a separate function and depend on the process
function PreFactorKPC(ptKinematic,scaleMu,proc1)
  type(DYpoint),intent(in)::ptKinematic
  integer,intent(in)::proc1
  real(dp),intent(in)::scaleMu
  real(dp)::PreFactorKPC

  SELECT CASE(proc1)
  case(0)
    PreFactorKPC=1d0
  CASE(1)
!     !4 pi aEm^2/3 /Nc/s
    PreFactorKPC=pi2x2/9*(alphaEM(scaleMu)**2)/ptKinematic%s*&
        HardCoefficientDY(scaleMu)*&
        hc2*1d9!from GeV to pb

  CASE(2)
    !4 pi aEm^2/3 /Nc/s
    ! the process=2 is for the xF-integration. It has extra weight J^{-1} with J= 2sqrt[(Q^2+q_T^2)/s] Cosh[y]
    PreFactorKPC=pi2x2/9*(alphaEM(scaleMu)**2)/ptKinematic%s*&
        HardCoefficientDY(scaleMu)*&
        hc2*1d9*&!from GeV to pb
        1._dp/(2._dp*ptKinematic%tau*coshY(ptKinematic)) !!!! J^{-1}

  CASE (5) !exclusive HIGGSboson production
    ! (2\pi) *pi Mh^2 as(mu)/36/s/vev^2 * H*cT^2
    ! (1.033)^2 is correction for mT mass in Ct at LO.
    ! * pi*Q2/2 (because of KPC ??? this factor is assumption)
    PreFactorKPC=(pi/36d0)*MH2*(As(c2_global*scaleMu)/VEVH)**2*ptKinematic%Q2/ptKinematic%s*&
        HardCoefficientHIGGS(scaleMu)*(EffCouplingHFF(scaleMu)**2)*1.0677023627519822d0*&
        hc2*1d9!from GeV to pb

  CASE (3,4)
    error stop ErrorString('Z/W in narrow-width approx. is not supported in KPC mode',moduleName)

  CASE DEFAULT
    error stop ErrorString('unknown process p='//numToStr(proc1)//' .Evaluation stop.',moduleName)
  END SELECT

end function PreFactorKPC

!!!! Fiducial cut for leptons.
!!!! this factor depends on the type of process and computed in the corresponding module
!!!! This function perform the preliminary selection depending on the process
function LeptonCutFactorKPC(ptKinematic,proc1, includeCuts_in,CutParam)
  type(DYpoint),intent(in)::ptKinematic
  logical,intent(in)::includeCuts_in
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,intent(in)::proc1
  real(dp)::LeptonCutFactorKPC

  if(includeCuts_in) then
  !!! here include cuts on the lepton tensor
  !!!!!!!!!!!!!!!
  !!!! the cases 2N correspond to part ~f1f1 in angular coefficient A_N
  !!!! the cases 3N correspond to part ~h1h1 in angular coefficient A_N
  SELECT CASE(proc1)
  CASE(1,2,4,5,6,101,102,103) !!! UU
      LeptonCutFactorKPC=CutFactor(qT_in=ptKinematic%qT,Q_in=ptKinematic%Q,y_in=ptKinematic%y,CutParameters=CutParam,Cut_Type=-1)

  CASE(20,30) !!! Angular coefficients A0
      LeptonCutFactorKPC=CutFactor(qT_in=ptKinematic%qT,Q_in=ptKinematic%Q,y_in=ptKinematic%y,CutParameters=CutParam,Cut_Type=0)

  CASE(21,31) !!! Angular coefficients  A1
      LeptonCutFactorKPC=CutFactor(qT_in=ptKinematic%qT,Q_in=ptKinematic%Q,y_in=ptKinematic%y,CutParameters=CutParam,Cut_Type=1)

  CASE(22,32) !!! Angular coefficients A2
      LeptonCutFactorKPC=CutFactor(qT_in=ptKinematic%qT,Q_in=ptKinematic%Q,y_in=ptKinematic%y,CutParameters=CutParam,Cut_Type=2)

  CASE(23) !!! Angular coefficients A3
      LeptonCutFactorKPC=CutFactor(qT_in=ptKinematic%qT,Q_in=ptKinematic%Q,y_in=ptKinematic%y,CutParameters=CutParam,Cut_Type=3)

  CASE(24) !!! Angular coefficients A4
      LeptonCutFactorKPC=CutFactor(qT_in=ptKinematic%qT,Q_in=ptKinematic%Q,y_in=ptKinematic%y,CutParameters=CutParam,Cut_Type=4)

  CASE(35) !!! Angular coefficients A5
      LeptonCutFactorKPC=CutFactor(qT_in=ptKinematic%qT,Q_in=ptKinematic%Q,y_in=ptKinematic%y,CutParameters=CutParam,Cut_Type=5)

  CASE(36) !!! Angular coefficients A6
      LeptonCutFactorKPC=CutFactor(qT_in=ptKinematic%qT,Q_in=ptKinematic%Q,y_in=ptKinematic%y,CutParameters=CutParam,Cut_Type=6)

  CASE DEFAULT
      !!! default is spherical cut-factor
      LeptonCutFactorKPC=CutFactor(qT_in=ptKinematic%qT,Q_in=ptKinematic%Q,y_in=ptKinematic%y,CutParameters=CutParam,Cut_Type=-1)

  END SELECT


  else
    !!! this is uncut lepton tensor
    LeptonCutFactorKPC=1._dp
  end if

end function LeptonCutFactorKPC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS CALCULATING CROSS-SECTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------UNINTEGRATED------------------------------------------------------------------

!!! Function which evaluates the cross-section of DY at single kinematical point
!!! Important note: it actually computed d\sigma/dY/dQ2/dqT2 (all Jacobians are defined in Prefactors)
function xSec_DY_1pt(ptKinematic,process,incCut,CutParam)
  real(dp):: xSec_DY_1pt,FF,LC
  real(dp)::x1,x2,scaleMu,scaleZeta
  type(DYpoint),intent(in)::ptKinematic
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process

  GlobalCounter=GlobalCounter+1
  CallCounter=CallCounter+1

  if(isConvergenceLost_inSubmodule()) then
    xSec_DY_1pt=1d9
    return
  end if

  !!!!!!!!!!!!!!!!!!!! COMPUTATION WITH KPC
  if(useKPC) then
    !!! setting values of X
    x1=ptKinematic%tau*ptKinematic%expy
    x2=ptKinematic%tau/ptKinematic%expy

    !!! setting values of scales
    !!! scales are always symmetric mu2=zeta
    !!! scaleZeta=ptKinematic%Q2  !! KPCs are defined only for zeta=mu2 so, no reason to define it.
    scaleMu=ptKinematic%Q

    FF=KPC_DYconv(ptKinematic%Q2,ptKinematic%qT,x1,x2,scaleMu*c2_global,process(2:4))
    LC=LeptonCutFactorKPC(ptKinematic,process(4),incCut,CutParam)
    xSec_DY_1pt=PreFactorKPC(ptKinematic,scaleMu,process(1))*FF*LC

  !!!!!!!!!!!!!!!!!!!! COMPUTATION WITHOUT KPC
  else
    if(exactX1X2) then
        x1=ptKinematic%tau*ptKinematic%expy
        x2=ptKinematic%tau/ptKinematic%expy
    else
        x1=ptKinematic%Q/sqrt(ptKinematic%s)*ptKinematic%expy
        x2=ptKinematic%Q/sqrt(ptKinematic%s)/ptKinematic%expy
    end if
    !!! setting values of scales
    !!! scales are always symmetric mu2=zeta
    if(exactScales) then
        !!!!!! setting the scale zeta of the factorization (zeta=Q2+qT^2)
        scaleZeta=ptKinematic%Q2+ptKinematic%qT*ptKinematic%qT
        scaleMu=sqrt(scaleZeta)
    else
        !!!!!! setting the scale zeta of the factorization Q
        scaleZeta=ptKinematic%Q2
        scaleMu=ptKinematic%Q
    end if


    !!!!! compute cross-section for each process
    FF=TMDF_F(ptKinematic%Q2,ptKinematic%qT,x1,x2,scaleMu*c2_global,scaleZeta,scaleZeta,process(2:4))
    LC=LeptonCutFactorLP(ptKinematic,process(4),incCut,CutParam)
    xSec_DY_1pt=PreFactor2(ptKinematic,scaleMu,process(1))*FF*LC
  end if
end function xSec_DY_1pt

end module TMDX_DY_point
