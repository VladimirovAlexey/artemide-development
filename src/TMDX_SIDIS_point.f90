!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.05
!
!    This is a part of TMDX_SIDIS complex of modules.
!    This module computes the cross-section of SIDIS reaction (and related sub-processes)
!    at a single kinematical point.
!
!    if you use this module please, quote 1912.06532
!
!    ver 3.05: created from v.3.04 (AV, 02.07.2026)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDX_SIDIS_point
use aTMDe_Numerics
use aTMDe_IO
use aTMDe_ptSpec
use TMDF
use TMDF_KPC
use QCDinput
use EWinput

implicit none
private

!Current version of module
character (len=14),parameter :: moduleName="TMDX-SIDIS-1pt"
character (len=5),parameter :: version="v3.05"
!Last appropriate version of constants-file
integer,parameter::inputver=39

real(dp) :: toleranceINT=0.0001d0
real(dp) :: toleranceGEN=0.0000001d0

integer::outputlevel=2
type(Warning_OBJ)::Warning_Handler

logical::started=.false.

real(dp)::hc2

!!other global parameters, which are defined upon initialization
integer:: orderH_global

!!!----------------------POWER CORRECTIONS---------------------------------
!!! The SIDIS cross-section is defined in laboratory frame, whereas the factorization theorem in Breit frame
!!! Thus, many variables have involved form, depending on masses of hadrons.
!!! However, the factorization theorem is strictly massless.
!!! Thus, if masses set non-zero, they are accounted in the kinematic factors (lepton tensor, Jacobians, etc),
!!! but not in the factorization variables (x,z,..)
!!!======= LP factorization
!! inclusion of power corrections to the definition
!! corrQT ~ qT/Q in cross-section
!! exactX1Z1 ~ x1,z1 as at LP factorization
!! exactZeta Use exact LP factorization scale =-2q^+q^-
logical:: corrQT,exactX1Z1, exactZeta !!!if true include
!!!======= KPC factorization
logical:: useKPC !!! use the KPC formula

real(dp)::c2_global

integer::GlobalCounter
integer::CallCounter


public::TMDX_SIDIS_1pt_Initialize,TMDX_SIDIS_1pt_ResetCounters,TMDX_SIDIS_1pt_SetScaleVariation
public::isConvergenceLost_inSubmodule,TMDX_SIDIS_1pt_IsInitialized
public::xSec_SIDIS_1pt




contains

function TMDX_SIDIS_1pt_IsInitialized()
    logical::TMDX_SIDIS_1pt_IsInitialized
    TMDX_SIDIS_1pt_IsInitialized=started
end function TMDX_SIDIS_1pt_IsInitialized

!! Initialization of the package
subroutine TMDX_SIDIS_1pt_Initialize(file,prefix)
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

    call MoveTO(51,'*0   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) FILEver
    if(FILEver<inputver) then
        CLOSE (51, STATUS='KEEP')
        error stop ErrorString('const-file version is too old. Update with artemide.setup.',moduleName)
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger

    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) hc2

    call MoveTO(51,'*10  ')
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
        CASE ("N3LO") !! same as NNNLO
            orderH_global=3
        CASE ("N4LO")
            orderH_global=4
        CASE DEFAULT
            if(outputLevel>0) write(*,*) WarningString('try to set unknown order. Switch to NLO.',moduleName)
            orderH_global=1
    END SELECT
    if(outputLevel>1) write(*,*) '    artemide.'//trim(moduleName)//': the used order is ',trim(orderMain)

    call MoveTO(51,'*p2  ')
    read(51,*) useKPC
    if(outputLevel>1 .and. useKPC) write(*,*) color('    artemide.'//trim(moduleName)//': using TMD factorization with KPC',c_cyan)
    if(outputLevel>1 .and. useKPC) write(*,*) color('                      Please, cite [2307.13054]',c_cyan)
    if(outputLevel>1 .and. .not.(useKPC)) &
      write(*,*) color('    artemide.'//trim(moduleName)//': using TMD factorization at LP',c_cyan)

    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) toleranceGEN
    call MoveTO(51,'*p2  ')
    read(51,*) toleranceINT


    !!------------------------------------LP FACTORIZATION

    call MoveTO(51,'*C   ')
    !! qT correction in kinematics
    call MoveTO(51,'*p1  ')
    read(51,*) corrQT
    if(outputLevel>2 .and. corrQT .and. (.not.useKPC)) &
            write(*,*) '    artemide.'//trim(moduleName)//': qT/Q corrections in kinematics are included.'
    !! qT correction in x1 z1
    call MoveTO(51,'*p4  ')
    read(51,*) exactX1Z1
    if(outputLevel>2 .and. exactX1Z1 .and. (.not.useKPC)) &
            write(*,*) '    artemide.'//trim(moduleName)//': Exact LP values for x1,z1 are included.'
    !!exact values for scales
    call MoveTO(51,'*p5  ')
    read(51,*) exactZeta
    if(outputLevel>2 .and. exactZeta .and. (.not.useKPC)) &
          write(*,*) '    artemide.'//trim(moduleName)//': Exact LP values of factorization scales variables are included.'

    !!------------------------------------KPC FACTORIZATION
    if(useKPC) then
      !!!! KPC must be computed with full expression for qT
      corrQT=.true.
      exactX1Z1=.true.
      exactZeta=.true.


    call MoveTO(51,'*D   ')
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
        if(outputLevel>1) write(*,*) '.. initializing TMDF_KPC (from ',moduleName,')'
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

    c2_global=1d0

    GlobalCounter=0
    CallCounter=0

    started=.true.
    write(*,*)  color('----- arTeMiDe.'//trim(moduleName)//' '//trim(version)//': .... initialized',c_green)
end subroutine TMDX_SIDIS_1pt_Initialize


subroutine TMDX_SIDIS_1pt_ResetCounters()
    if(outputlevel>2) then
        write(*,'(A,ES12.3)') 'TMDX SIDIS statistics   total calls of xSec (complete session) :  ',Real(GlobalCounter)
        write(*,'(A,ES12.3)') '                        total calls of xSec (current NP-setup) :  ',Real(CallCounter)
    end if
    CallCounter=0
    call Warning_Handler%Reset()
end subroutine TMDX_SIDIS_1pt_ResetCounters


!!!!Call this after TMD initialization but before NP, and X parameters
subroutine TMDX_SIDIS_1pt_SetScaleVariation(c2_in)
real(dp),intent(in)::c2_in

if(outputLevel>1) write(*,*) trim(moduleName)//': scale variation constant c2 reset:',c2_in

if(c2_in<0.1d0 .or. c2_in>10.d0) then
    if(outputLevel>0) write(*,*) WarningString('variation in c2 is enormous. c2 is set to 2',moduleName)
    c2_global=2d0
else
    c2_global=c2_in
end if

end subroutine TMDX_SIDIS_1pt_SetScaleVariation

!!!!! wrapper, to check the convergence in proper module.
function isConvergenceLost_inSubmodule()
logical::isConvergenceLost_inSubmodule
if(useKPC) then
  isConvergenceLost_inSubmodule=TMDF_KPC_IsconvergenceLost()
else
  isConvergenceLost_inSubmodule=TMDF_IsconvergenceLost()
end if
end function isConvergenceLost_inSubmodule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS FOR PREFACTORS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! for easier read coeff-functions are split into separate file
INCLUDE '/Code/TMDX/SIDIScoeff-func.f90'
  
!!!!! Prefactor2 is (universal part) x H
!!!!! In addition to usual p and processes, it receives the following variables (to be computed only once, to save time)
!!!!! x1, z1 = LP variables x1 and z1
!!!!! qT = photon's transverse momentum
!!!!! eps = SIDIS variables varepsilon
!!!!! scaleMU = factorization scale MU
function PreFactor2(p,process,x1,z1,qT,eps,scaleMu)
type(SIDISpoint),intent(in)::p
integer,dimension(1:4),intent(in)::process
real(dp),intent(in)::x1,z1,qT,eps,scaleMu
real(dp)::PreFactor2,uniPart,fac1,jacob

!!!! Note, that these factors have extra factor 2, which is compensated by the definition of Fourier transform F with 1/2 (in TMDF)

jacob=pJacobian(p)

!! fac1 is the prefactor for the unpolarized expression, which appears at LP due to (-g_T*W)
!!! this is 1+qT^2/Q^2(e-gamma^2/2)/(1+gamma^2)
if(corrQT) then
    fac1=1d0+(qT**2/p%Q2)*(eps-0.5d0*p%gamma2)/(1+p%gamma2)
else
    fac1=1._dp
end if

SELECT CASE(process(1))
  case(0)
      PreFactor2=1._dp
  CASE(1)
      !!! uniPart is the prefactor for the cross-section
      !2 pi aEm^2/Q^4 y^2/(1-epsilon)/J
      uniPart=pix2*alphaEM(scaleMu)**2/(p%Q2**2)*(p%y**2)/((1d0-eps)*jacob)

      PreFactor2=uniPart*fac1*HardCoefficientSIDIS(scaleMu)*&
      hc2*1d9!from GeV to pbarn

  CASE(2)
      !!!! same as for case(1) but with factor Q2/y
      uniPart=pix2*alphaEM(scaleMu)**2/(p%Q2)*(p%y)/((1d0-eps)*jacob)

      PreFactor2=uniPart*fac1*HardCoefficientSIDIS(scaleMu)*&
      hc2*1d9!from GeV to pbarn

  CASE(3)
      !!!! same as for case(1) but with factor x/y
      uniPart=pix2*alphaEM(scaleMu)**2/(p%Q2**2)*(p%y*p%x)/((1d0-eps)*jacob)

      PreFactor2=uniPart*fac1*HardCoefficientSIDIS(scaleMu)*&
      hc2*1d9!from GeV to pbarn

  CASE DEFAULT
      PreFactor2=0._dp
      error stop ErrorString(' unknown process p1='//numToStr(process(1)),moduleName)
END SELECT
end function PreFactor2


!!!!! Prefactor for KPC-case is (universal part) x H
!!!!! In addition to usual p and processes, it receives the following variables (to be computed only once, to save time)
!!!!! x1, z1 = LP variables x1 and z1
!!!!! qT = photon's transverse momentum
!!!!! eps = SIDIS variables varepsilon
!!!!! scaleMU = factorization scale MU
function PreFactorKPC(p,process,x1,z1,qT,eps,scaleMu)
type(SIDISpoint),intent(in)::p
integer,dimension(1:4),intent(in)::process
real(dp),intent(in)::x1,z1,qT,eps,scaleMu
real(dp)::PreFactorKPC,uniPart,jacob

!!!-----------------------------------------------------------------
!!!------  Prefactors KPC case (TO BE CHECKED!!) [there is a suspicious factor 1/z]
!!!-----------------------------------------------------------------

jacob=pJacobian(p)

SELECT CASE(process(1))
  case(0)
      PreFactorKPC=1._dp
  CASE(1)
      !!! uniPart is the prefactor for the cross-section
      !pi^2 aEm^2/Q^2 y^2/(1-epsilon)/J
      uniPart=pi2*alphaEM(scaleMu)**2/(p%Q2)*(p%y)**2/((1d0-eps)*jacob)

      PreFactorKPC=uniPart*HardCoefficientSIDIS(scaleMu)*&
      hc2*1d9!from GeV to pbarn

  CASE(2)
      !!!! same as for case(1) but with factor Q2/y

      uniPart=pi2*alphaEM(scaleMu)**2*(p%y)/((1d0-eps)*jacob)

      PreFactorKPC=uniPart*HardCoefficientSIDIS(scaleMu)*&
      hc2*1d9!from GeV to pbarn

  CASE(3)
      !!!! same as for case(1) but with factor x/y

      !!! uniPart is the prefactor for the cross-section
      uniPart=pi2*alphaEM(scaleMu)**2/(p%Q2)*(p%y*p%x)/((1d0-eps)*jacob)

      PreFactorKPC=uniPart*HardCoefficientSIDIS(scaleMu)*&
      hc2*1d9!from GeV to pbarn

  CASE DEFAULT
      PreFactorKPC=0._dp
      error stop ErrorString(' unknown process p1='//numToStr(process(1)),moduleName)
END SELECT
end function PreFactorKPC
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS CALCULATING CROSS-SECTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------UNINTEGRATED------------------------------------------------------------------

!!! Evaluates the SIDIS cross-section at a single kinematic point
function xSec_SIDIS_1pt(p,process)
type(SIDISpoint),intent(in)::p
integer,dimension(1:4),intent(in)::process

real(dp):: xSec_SIDIS_1pt,FF
real(dp)::x1,z1,qT,eps,scaleMu,scaleZeta


GlobalCounter=GlobalCounter+1
CallCounter=CallCounter+1

if(isConvergenceLost_inSubmodule()) then
  xSec_SIDIS_1pt=1d9
  return
end if

eps=eps_SIDIS(p)

!!!! KPC formula
if(useKPC) then
  call ComputeX1Z1qT(p,x1,z1,qT)
  scaleMu=p%Q
  scaleZeta=p%Q2

  FF=KPC_SIDISconv(p%Q2,qT,x1,z1,eps,scaleMu*c2_global,process(2:4))
  xSec_SIDIS_1pt=PreFactorKPC(p,process,x1,z1,qT,eps,scaleMu)*FF

else !!!! LP formula
  if(exactX1Z1) then
    !!!! Note, this subroutine, automatically provide exact LP values
    call ComputeX1Z1qT(p,x1,z1,qT)
  else
    !!! use simplified formulas
    qT=qT_SIDIS(p)
    x1=p%x
    z1=p%z
  end if

  !!! setting values of scales
  !!!! If exact scales, zeta*zeta=Q^2-qT^2.
  if(exactZeta) then
    scaleZeta=p%Q2-qT**2
    scaleMu=sqrt(scaleZeta)
  else
    scaleMu=p%Q
    scaleZeta=p%Q2
  end if


  FF=TMDF_F(p%Q2,qT,x1,z1,scaleMu*c2_global,scaleZeta,scaleZeta,process(2:4))
  xSec_SIDIS_1pt=PreFactor2(p,process,x1,z1,qT,eps,scaleMu)*FF

end if

end function xSec_SIDIS_1pt

end module TMDX_SIDIS_point
