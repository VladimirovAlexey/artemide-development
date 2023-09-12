!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.00
!
!    Evaluation of the TMD cross-section for DY-like cross-sections
!
!    if you use this module please, quote 1706.01473
!
!    ver 3.0: created from v.2.06 (AV, 05.09.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDX_DY
use aTMDe_Numerics
use IO_functions
use IntegrationRoutines
use TMDF
use LeptonCutsDY
use QCDinput
use EWinput

implicit none
private

!Current version of module
character (len=7),parameter :: moduleName="TMDX-DY"
character (len=5),parameter :: version="v3.00"
!Last appropriate verion of constants-file
integer,parameter::inputver=30

real(dp) :: toleranceINT=0.0001d0
real(dp) :: toleranceGEN=0.0001d0

integer::outputlevel
integer::messageTrigger

!!!!
!!!! in the module the kinematic is stored in the varibles "kinematic" real(dp),dimension(1:6)
!!!! which is (qT,s,Q,Q^2,x0,y,exp[y])
!!!! where x0=sqrt[(Q^2+q_T^2)/s]   (if exactX1X2) or x0=Q^2/s (otherwise)
!!!!

!!other global parameters see SetXParameters
integer:: orderH_global
logical::usePIresum
integer:: exactX1X2    !!!=1 if exact x's=true, =0 otherwise
integer:: exactScales  !!!=1 if exact hard scales = true, =0 otherwise

!!! number of sections for PT-integral by default
integer::NumPTdefault=4
!!! Maximum size of Q-bin. Larger bins are desected
real::maxQbinSize=30.

real(dp)::c2_global



integer::GlobalCounter
integer::CallCounter
integer::messageCounter

real(dp)::hc2

logical::started=.false.

public::TMDX_DY_Initialize,TMDX_DY_SetScaleVariation,TMDX_DY_ShowStatistic,TMDX_DY_ResetCounters,TMDX_DY_IsInitialized
public::  xSec_DY,xSec_DY_List,xSec_DY_List_BINLESS

interface xSec_DY
  module procedure MainInterface_AsAAAloo
end interface
 
contains
  
function TMDX_DY_IsInitialized()
  logical::TMDX_DY_IsInitialized
  TMDX_DY_IsInitialized=started
end function TMDX_DY_IsInitialized

!! Initialization of the package
subroutine TMDX_DY_Initialize(file,prefix)
  character(len=*),intent(in)::file
  character(len=*),intent(in),optional::prefix
  character(len=300)::path
  logical::initRequired,dummyLogical
  character(len=8)::orderMain
  integer::i,FILEver
  !$ integer:: omp_get_thread_num

  if(started) return

  if(present(prefix)) then
    path=trim(adjustl(prefix))//trim(adjustr(file))
  else
    path=trim(adjustr(file))
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
    stop
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

  !$    if(outputLevel>1) write(*,*) '    artemide.TMDX_DY: parallel evaluation of cross-sections is to be used'
  !$    call MoveTO(51,'*C   ')
  !$    call MoveTO(51,'*p1  ')
  !$    read(51,*) i
  !$    call OMP_set_num_threads(i)
  !$    if(outputLevel>1) write(*,*) '    artemide.TMDX_DY: number of threads for parallel evaluation is set to ', i

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

  !!exact values of x1x2
  call MoveTO(51,'*p2   ')
  read(51,*) dummyLogical
  if(dummyLogical) then
    exactX1X2=1
  else
    exactX1X2=0
  end if

  if(outputLevel>2 .and. dummyLogical) write(*,*) '    artemide.TMDX_DY: qT/Q corrections for x1 and x2 variables are included.'
  !! pi2 resummation
  call MoveTO(51,'*p3   ')
  read(51,*) usePIresum
  if(outputLevel>2 .and. usePIresum) write(*,*) '    artemide.TMDX_DY: pi-resummation in coef.function included.'
  !!exact values for scales
  call MoveTO(51,'*p4   ')
  read(51,*) dummyLogical
  if(dummyLogical) then
    exactScales=1
  else
    exactScales=0
  end if

  if(outputLevel>2 .and. dummyLogical) write(*,*) '    artemide.TMDX_DY: qT/Q correction for scales variables are included.'

  call MoveTO(51,'*B   ')
  call MoveTO(51,'*p1  ')
  read(51,*) toleranceGEN
  call MoveTO(51,'*p2  ')
  read(51,*) toleranceINT
  call MoveTO(51,'*p3  ')
  read(51,*) NumPTdefault
  call MoveTO(51,'*p4  ')
  read(51,*) maxQbinSize


  !$     if(outputLevel>2) write(*,*) '------TEST OF PARALLEL PROCESSING ----------'
  !$OMP PARALLEL
  !$     if(outputLevel>2) write(*,*) '   artemide.TMDX_DY:thread num ',  omp_get_thread_num(), ' ready.'
  !$OMP END PARALLEL
  CLOSE (51, STATUS='KEEP')

  if(.not.EWinput_IsInitialized()) then
    if(outputLevel>1) write(*,*) '.. initializing EWinput (from ',moduleName,')'
    if(present(prefix)) then
      call EWinput_Initialize(file,prefix)
    else
      call EWinput_Initialize(file)
    end if
  end if

  if(.not.TMDF_IsInitialized()) then
    if(outputLevel>1) write(*,*) '.. initializing TMDF (from ',moduleName,')'
    if(present(prefix)) then
      call TMDF_Initialize(file,prefix)
    else
      call TMDF_Initialize(file)
    end if
  end if

  c2_global=1d0

  GlobalCounter=0
  CallCounter=0
  messageCounter=0

  started=.true.

  write(*,*)  color('----- arTeMiDe.TMD_DY '//trim(version)//': .... initialized',c_green)
end subroutine TMDX_DY_Initialize

!!! routine that reset counter to 0
subroutine TMDX_DY_ResetCounters()
  if(outputlevel>2) call TMDX_DY_ShowStatistic
  GlobalCounter=0
  CallCounter=0
  messageCounter=0
end subroutine TMDX_DY_ResetCounters

!!! Show current state
subroutine TMDX_DY_ShowStatistic()
    write(*,'(A,ES12.3)') 'TMDX DY statistics      total calls of point xSec  :  ',Real(GlobalCounter)
    write(*,'(A,ES12.3)') '                              total calls of xSecF :  ',Real(CallCounter)
    write(*,'(A,F12.3)')  '                                         avarage M :  ',Real(GlobalCounter)/Real(CallCounter)
end subroutine TMDX_DY_ShowStatistic
  
  
  
!!! Vary the scale c2
subroutine TMDX_DY_SetScaleVariation(c2_in)
  real(dp),intent(in)::c2_in

  if(outputLevel>1) write(*,*) 'TMDX_DY: c2 scale reset:',c2_in

  if(c2_in<0.1d0 .or. c2_in>10.d0) then
    if(outputLevel>0) write(*,*) WarningString('variation in c2 is enourmous. c2 is set to 2',moduleName)
    c2_global=2d0
  else
    c2_global=c2_in
  end if

end subroutine TMDX_DY_SetScaleVariation


  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PROCESS DEFINITION:
!!! (i, n, h1, h2) -- integers
!!! i = the integration type and global prefactor
!!! n = the integrand combination
!!! h1,h2 = the hadron types. h>0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! this function returns .true. is the process is symmetric with repsect y->-y.
!!! this is used to speed up the computation of symmetric cross-section in the symmetric regions.
function IsySymmetric(process)
  integer,intent(in),dimension(1:4)::process
  logical:: IsySymmetric

  !!! hadrons should be same and it is not xF-integral
  if(process(3)/=process(4) .or. process(1)==2) then
    IsySymmetric=.false.
    return
  end if

  !!! only few processes are symetric these are DY, W,Z, H.
  if(process(2)<=22) then
    IsySymmetric=.true.
    return
  end if
  IsySymmetric=.false.

end function IsySymmetric


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS FOR OPERATION WITH KINEMATICS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
!!! function makes kinematic array from the given set of qT,s,Q,y
!!! array has 7 often appearing entries
pure function kinematicArray(qT,s,Q,y)
real(dp),dimension(1:7)::kinematicArray
real(dp),intent(in)::qT,s,Q,y

  kinematicArray=(/qT,s,Q,Q**2,sqrt((Q**2+exactX1X2*qT**2)/s),y,exp(y)/)

end function kinematicArray
  
!!!intrinsic change the value of Q within kinematic array var
subroutine SetQ(Q,var)
  real(dp),dimension(1:7)::var
  real(dp),intent(in)::Q

  var(3)=Q
  var(4)=Q**2
  var(5)=sqrt((var(4)+exactX1X2*var(1)**2)/var(2))

end subroutine SetQ
  
!!!intrinsic change the value of y within kinematic array var
subroutine SetY(y,var)
  real(dp),dimension(1:7)::var
  real(dp),intent(in)::y

  var(6)=y
  var(7)=exp(y)

end subroutine SetY
  
!!!intrinsic change the value of qT within kinematic array var
subroutine SetQT(qT_in,var)
  real(dp),dimension(1:7)::var
  real(dp),intent(in)::qT_in

  var(1)=qT_in
  var(5)=sqrt((var(4)+exactX1X2*var(1)**2)/var(2))

end subroutine SetQT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS FOR PREFACTORS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! for easier read coeff-functions are split into separate file
INCLUDE 'DYcoeff-func.f90'
  
!!!!! Prefactor 2 is (universal part) x (cuts) x H
function PreFactor2(kin,process, includeCuts_in,CutParam)
  real(dp),dimension(1:7),intent(in)::kin
  logical,intent(in)::includeCuts_in
  real(dp)::PreFactor2,cutPrefactor,uniPart,scaleMu
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process

  !!!!! cut part
  if(includeCuts_in) then
      !!! here include cuts onf lepton tensor
      cutPrefactor=CutFactor4(qT=kin(1),Q_in=kin(3),y_in=kin(6),CutParameters=CutParam)
  else
  !!! this is uncut lepton tensor
      cutPrefactor=(1+0.5_dp*(kin(1)/kin(3))**2)
  end if


  !!!! universal part
  scaleMu=sqrt(kin(4)+exactScales*kin(1)**2)

SELECT CASE(process(1))
  case(-10221191)
  uniPart=1d0
  cutPrefactor=1d0
  CASE(1)
  !4 pi aEm^2/3 /Nc/Q^2/s
  uniPart=pix4/9d0*(alphaEM(scaleMu)**2)/(kin(2)*kin(4))*&
      HardCoefficientDY(scaleMu)*&
      hc2*1d9!from GeV to pb

  CASE(2)
  !4 pi aEm^2/3 /Nc/Q^2/s
  ! the process=2 is for the xF-integration. It has extra weigth 2sqrt[(Q^2+q_T^2)/s] Cosh[y]
  uniPart=pix4/9d0*(alphaEM(scaleMu)**2)/(kin(2)*kin(4))*&
      HardCoefficientDY(scaleMu)*&
      hc2*1d9*&!from GeV to pb
      2._dp*kin(5)*cosh(kin(6))

  CASE (3) !Zboson in the narrow-width approximation
  !4 pi^2 aem/Nc/s Br(z->ee+mumu)
  uniPart=pi2x4/3d0*alphaEM(scaleMu)/kin(2)*&
      HardCoefficientDY(scaleMu)*&
      hc2*1d9*&!from GeV to pb
      0.03645d0!Br from PDG, ee+mumu

  CASE (4) !Wboson in the narrow-width approximation
  !4 pi^2 aem/Nc/s Br(z->ee+mumu)
  uniPart=pi2x4/3d0**alphaEM(scaleMu)/kin(2)*&
      HardCoefficientDY(scaleMu)*&
      hc2*1d9*&!from GeV to pb
      0.1086d0!Br from PDG, ee+mumu

  CASE (5) !exclusive HIGGSboson production
  ! (2\pi) *pi Mh^2 as(mu)/36/s/vev^2 * H*cT^2
  ! (1.033)^2 is correction for mT mass in Ct at LO.
  uniPart=(1d0/18d0)*MH2*(As(c2_global*scaleMu)/VEVH)**2/kin(2)*&
      HardCoefficientHIGGS(scaleMu)*(EffCouplingHFF(scaleMu)**2)*1.0677023627519822d0*&
      hc2*1d9!from GeV to pb

  cutPrefactor=1d0 !!! cut-prefactor is different in this case!
  CASE DEFAULT
    write(*,*) 'ERROR: arTeMiDe.TMDX_DY: unknown process p2=',process(1),' .Evaluation stop.'
    stop
END SELECT

PreFactor2=uniPart*cutPrefactor

end function PreFactor2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS CALCULATING CROSS-SECTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!---------------------------------UNINTEGRATED------------------------------------------------------------------

!!! this is the base function which evaluate xSec at single qt (without lists)
!!! this is extended (and default) version of xSec, which include all parameters
function xSec(var,process,incCut,CutParam)
  real(dp):: xSec,FF
  real(dp)::x1,x2,scaleMu,scaleZeta
  real(dp),dimension(1:7),intent(in)::var
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process
  GlobalCounter=GlobalCounter+1

  if(TMDF_IsconvergenceLost()) then
    xSec=1d9
    return
  end if

  !!! setting values of X
  x1=var(5)*var(7)
  x2=var(5)/var(7)


  !!! setting values of scales
  scaleZeta=var(4)+exactScales*var(1)**2  !! zeta=Q2+qT^2
  scaleMu=sqrt(scaleZeta)

  FF=TMDF_F(var(4),var(1),x1,x2,scaleMu*c2_global,scaleZeta,scaleZeta,process(2))

  xSec=PreFactor2(var,process,incCut,CutParam)*FF

  !write(*,*) "{",var(4),",",x1,"},"!,z1
end function xSec
  
  !---------------------------------INTEGRATED over Y---------------------------------------------------------------
  
!!! function determines the best value of PT-sections from PT-bin size, and Q
!!! it is determined by formula Q/PT< val/ (2 k) => def+2K
function NumPT_auto(dPT,Q)
  real,parameter::val=40.
  real(dp)::dPT,Q,rat
  integer::i,NumPT_auto
  rat=Q/dPT

  if(rat>40.) then
      NumPT_auto=NumPTdefault
      return
  else
      do i=1,25
          if(rat>(40d0/2d0/i)) then
              NumPT_auto=NumPTdefault+2*i
              return
          end if
      end do
  end if
  if(outputlevel>1) then
  write(*,*) WarningString('Fail to automatically determine number of Pt-section for a bin.',moduleName)
  write(*,*) '  ... Possibly Pt-bin is too large', dPT
  end if
  NumPT_auto=NumPTdefault+12

end function NumPT_auto

!!! Computes y-variable from x_F
pure function yFromXF(xF,var)
real(dp),dimension(1:7),intent(in)::var
real(dp), intent(in):: xF
real(dp):: yFromXF
  yFromXF=asinh(xF/2._dp/var(5))
end function yFromXF

!!!!!!!! This function computes the integration over (y,Q,pT)
!!!!!!!! it check the boundaries and split the integral and send to _IntegrateALL_0 function
function Xsec_IntegrateALL(process,incCut,CutParam,s_in,ymin_in,ymax_in,Qmin_in,Qmax_in,qtMin_in,qtMax_in)
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  real(dp) :: Xsec_IntegrateALL
  real(dp),intent(in):: s_in,ymin_in,ymax_in,Qmin_in,Qmax_in,qtMin_in,qtMax_in
  real(dp) :: s,Qmax,Qmin,qTmin,qTmax
  integer,dimension(1:4),intent(in)::process
  real(dp),dimension(1:7)::var
  integer::i,numSecPT,numSecQ

  if(TMDF_IsconvergenceLost()) then
    Xsec_IntegrateALL=1d9
    return
  end if

  !!!------------------------- checking Q----------
  if(Qmin_in<0.9d0) then
    call Warning_Raise('Attempt to compute xSec with Q<0.9.',messageCounter,messageTrigger,moduleName)
    write(*,*) "Qmin =",Qmin_in," (Qmin set to 1.GeV)"
    Qmin=1._dp
  else
    Qmin=Qmin_in
  end if

  if(Qmax_in<Qmin) then
    call Warning_Raise('Attempt to compute xSec with Qmax<Qmin. RESULT 0',messageCounter,messageTrigger,moduleName)
    Xsec_IntegrateALL=0._dp
    return
  end if
  Qmax=Qmax_in

  !!!------------------------- checking PT----------
  if(qTmin_in<0.0d0) then
    call Warning_Raise('Attempt to compute xSec with qT<0.',messageCounter,messageTrigger,moduleName)
    write(*,*) "qTmin =",qTmin_in," (qTmin set to 0.GeV)"
    qTmin=1._dp
  else
    qTmin=qTmin_in
  end if

  if(qTmax_in<qTmin) then
    call Warning_Raise('Attempt to compute xSec with qTmax<qTmin. RESULT 0',messageCounter,messageTrigger,moduleName)
    Xsec_IntegrateALL=0._dp
    return
  end if
  qTmax=qTmax_in

  !!!------------------------- checking S----------
  if(s_in<0.9d0) then
    call Warning_Raise('Attempt to compute xSec with s<0.9.',messageCounter,messageTrigger,moduleName)
    write(*,*) "s =",s_in," (s set to Qmin)"
    s=Qmin
  else
    s=s_in
  end if

  var=kinematicArray(qTmax,s,Qmax,yMax_in)


  !!! check how many PT-bins is here
  numSecPT=NumPT_auto(qtMax-qtMin,Qmin)

  Xsec_IntegrateALL=Xsec_IntegrateALL_0(var,process,incCut,CutParam,yMin_in,yMax_in,Qmin,Qmax,qtMin,qtMax,numSecPT)

end function Xsec_IntegrateALL

!!!!! before this function all parameers must be checked
!!!!! This is triple-iterative function
!!!!! it has 3 included definition of functions....
!!!!! The y is SA, the Q is SA, the PT is SN
function Xsec_IntegrateALL_0(var,process,incCut,CutParam,ymin_in,ymax_in,Qmin_in,Qmax_in,qtMin_in,qtMax_in,numPT)
  real(dp),dimension(1:7) :: var
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  real(dp),intent(in):: ymin_in,ymax_in,Qmin_in,Qmax_in,qtMin_in,qtMax_in
  integer,dimension(1:4),intent(in)::process
  integer,intent(in)::numPT
  real(dp)::Xsec_IntegrateALL_0,qtCUR


  real(dp)::delta,inter,qT_cur
  integer::i

  !!!!!!!!!!!!!!!!!!!fixed number Simpsons
    delta=(qTMax_in-qTMin_in)/numPT
    call SetQT(qTMin_in,var)
    inter=2*qTMin_in*Xsec_Qint_Yint(var,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in)!!!! first term

    !!!! even terms
    do i=1,numPT-1,2
        qtCUR=qMin_in+i*delta
        call SetQT(qtCUR,var)
        inter=inter+8*qtCUR*Xsec_Qint_Yint(var,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in)
    end do
    !!!! odd term
    do i=2,numPT-2,2
        qtCUR=qMin_in+i*delta
        call SetQT(qtCUR,var)
        inter=inter+4*qtCUR*Xsec_Qint_Yint(var,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in)
    end do

    call SetQT(qMax_in,var)
    inter=inter+2*qMax_in*Xsec_Qint_Yint(var,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in)!!!! last term

    Xsec_IntegrateALL_0=delta/3._dp*inter
end function Xsec_IntegrateALL_0

!!!
function Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
  real(dp),dimension(1:7) :: var
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  real(dp) :: Xsec_Yint
  real(dp) :: ymin, ymax,ymin_in,ymax_in
  real(dp) :: ymin_Check,ymax_Check
  integer,dimension(1:4),intent(in)::process

  if(TMDF_IsconvergenceLost()) then
    Xsec_Yint=1d9
    return
  end if

  !!! evaluate correspnding y's
  !!! in the case process=2 the integral is over xF
  if(process(1)==2) then
    ymin=yFromXF(ymin_in,var)
    ymax=yFromXF(ymax_in,var)
  else
    ymin=ymin_in
    ymax=ymax_in
  end if

  ymin_Check=log(var(5))+0.000000001d0
  ymax_Check=-log(var(5))-0.000000001d0

  if(IsySymmetric(process) .and. (ABS(ymax+ymin)<toleranceGEN)) then!!! symetric integral
  if(ymax > ymax_check) then
      ymax=ymax_Check
  end if!!!!! else case: automatically taken into account

  Xsec_Yint=2d0*integralOverYpoint_S(var,process,incCut,CutParam,0d0,ymax)!!! 2 since symmetric

  else !!!non-symmetric integral!!!!!!!!
    if(ymax<ymin_check .or. ymin>ymax_check) then !!! the case then y is outside physicsl region
      Xsec_Yint=0d0
    else
      if(ymax > ymax_check) then
        ymax=yMax_check
      end if!!!!! else case: automatically taken into account
      if(ymin < ymin_check) then
        ymin=ymin_check
      end if!!!!! else case: automatically taken into account

      Xsec_Yint=integralOverYpoint_S(var,process,incCut,CutParam,ymin,ymax)

    end if
  end if
end function Xsec_Yint
!--------------Simpsons--------------------
!!!! parameter valueMax remembers the approximate value of integral to weight the tolerance.
!!!! evaluation is done by adaptive simpson
!!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
!!!! Thus minimal number of points =9
function integralOverYpoint_S(var,process,incCut,CutParam,yMin_in,yMax_in)
  real(dp),dimension(1:7)::var
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process
  real(dp) ::integralOverYpoint_S
  real(dp) :: X1,X2,X3,X4,X5
  real(dp) :: y2,y3,y4,deltay
  real(dp) :: yMin_in,yMax_in
  real(dp)::valueMax

  deltay=(yMax_in-yMin_in)/4._dp
  y2=yMin_in+deltay
  y3=yMin_in+deltay*2
  y4=yMax_in-deltay

  call SetY(yMin_in,var)
  X1= xSec(var,process,incCut,CutParam)
  call SetY(y2,var)
  X2= xSec(var,process,incCut,CutParam)
  call SetY(y3,var)
  X3= xSec(var,process,incCut,CutParam)
  call SetY(y4,var)
  X4= xSec(var,process,incCut,CutParam)
  call SetY(yMax_in,var)
  X5= xSec(var,process,incCut,CutParam)

  !!approximate integral value
  valueMax=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0

  integralOverYpoint_S=IntegralOverYpoint_S_Rec(var,process,incCut,CutParam,yMin_in,y3,X1,X2,X3,valueMax)+&
    IntegralOverYpoint_S_Rec(var,process,incCut,CutParam,y3,yMax_in,X3,X4,X5,valueMax)
end function integralOverYpoint_S

!!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
recursive function integralOverYpoint_S_Rec(var,process,incCut,CutParam,yMin_in,yMax_in,X1,X3,X5,valueMax) result(interX)
  real(dp),dimension(1:7) ::var
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process
  real(dp) :: interX,X1,X2,X3,X4,X5
  real(dp) :: valueAB,valueACB
  real(dp) :: yMin_in,yMax_in,y2,y3,y4,deltay
  real(dp),intent(in)::valueMax

  deltay=(yMax_in-yMin_in)/4._dp
  y2=yMin_in+deltay
  y3=yMin_in+deltay*2
  y4=yMax_in-deltay


  call SetY(y2,var)
  X2= xSec(var,process,incCut,CutParam)

  call SetY(y4,var)
  X4= xSec(var,process,incCut,CutParam)


  valueAB=deltay*(X1+4d0*X3+X5)/6d0
  valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0

  If(ABS((valueACB-valueAB)/valueMax)>toleranceINT) then
  interX=integralOverYpoint_S_Rec(var,process,incCut,CutParam,yMin_in,y3,X1,X2,X3,valueMax)&
    +integralOverYpoint_S_Rec(var,process,incCut,CutParam,y3,yMax_in,X3,X4,X5,valueMax)
  else
  interX=valueACB
  end if

end function integralOverYpoint_S_Rec

!---------------------------------INTEGRATED over Y over Q---------------------------------------------------------------
!!!! No need for check over Y they take a place within y-integration

!!!! to integrate over Q I use adaptive Simpson. (defined below)
!!!! before the integration over Q I check the size of Q-bin,
!!!! if Q-bin is large I split desect the integration range to smaller
function Xsec_Qint_Yint(var,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in)
  real(dp),dimension(1:7)::var
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process
  real(dp),intent(in) :: yMin_in,yMax_in,QMin_in,QMax_in
  real(dp):: Xsec_Qint_Yint
  integer::numSec,i
  real(dp)::dQ

  if(TMDF_IsconvergenceLost()) then
    Xsec_Qint_Yint=1d9
    return
  end if

  !!! check how many maxQbins is inside the integration range (+1)
  numSec=INT((Qmax_in-Qmin_in)/maxQbinSize)+1

  if(numSec==1) then
    !!! if the bin is smaller than maxQbinSize, integrate as is
    Xsec_Qint_Yint=Xsec_Qint_Yint_in(var,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in)
  else
    !!! else divide to smaler bins and sum the integrals
    dQ=(Qmax_in-Qmin_in)/numSec !!! size of new bins

    Xsec_Qint_Yint=0d0
    do i=0,numSec-1
      Xsec_Qint_Yint=Xsec_Qint_Yint + &
          Xsec_Qint_Yint_in(var,process,incCut,CutParam,Qmin_in+i*dQ,Qmin_in+(i+1)*dQ,ymin_in,ymax_in)
    end do
  end if
end function Xsec_Qint_Yint

!--------------Simpsons--------------------
!!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
!!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
!!!! Thus minimal number of points =9
!!!! taking into account minimum calls of y-integral we have  =81 points
function Xsec_Qint_Yint_in(var,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in)
real(dp),dimension(1:7)::var
logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process
  real(dp),intent(in) :: yMin_in,yMax_in,QMin_in,QMax_in
  real(dp):: Xsec_Qint_Yint_in
  real(dp) :: X1,X2,X3,X4,X5
  real(dp)::valueMax,Q2,Q3,Q4,deltaQ

  deltaQ=(QMax_in-QMin_in)/4._dp
  Q2=QMin_in+deltaQ
  Q3=QMin_in+deltaQ*2
  Q4=QMax_in-deltaQ

  call SetQ(QMin_in,var)
  X1=2*QMin_in*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)

  call SetQ(Q2,var)
  X2=2*Q2*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)

  call SetQ(Q3,var)
  X3=2*Q3*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)

  call SetQ(Q4,var)
  X4=2*Q4*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)

  call SetQ(QMax_in,var)
  X5=2*QMax_in*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)

    !!approximate integral value
  valueMax=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0

  Xsec_Qint_Yint_in=IntegralOverQYpoint_S_Rec(var,process,incCut,CutParam,QMin_in,Q3,yMin_in,yMax_in,X1,X2,X3,valueMax)+&
  IntegralOverQYpoint_S_Rec(var,process,incCut,CutParam,Q3,QMax_in,yMin_in,yMax_in,X3,X4,X5,valueMax)
end function Xsec_Qint_Yint_in

!!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
recursive function integralOverQYpoint_S_Rec(var,process,incCut,CutParam,&
                QMin_in,QMax_in,yMin_in,yMax_in,X1,X3,X5,valueMax) result(interX)
  real(dp),dimension(1:7)::var
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process
  real(dp) :: interX,X1,X2,X3,X4,X5
  real(dp) :: valueAB,valueACB
  real(dp) :: yMin_in,yMax_in,QMin_in,QMax_in,Q2,Q3,Q4,deltaQ
  real(dp),intent(in)::valueMax

  deltaQ=(QMax_in-QMin_in)/4._dp
  Q2=QMin_in+deltaQ
  Q3=QMin_in+deltaQ*2
  Q4=QMax_in-deltaQ

  call SetQ(Q2,var)
  X2=2*Q2*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)

  call SetQ(Q4,var)
  X4=2*Q4*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)

  valueAB=deltaQ*(X1+4d0*X3+X5)/6d0
  valueACB=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0

  If(ABS((valueACB-valueAB)/valueMax)>toleranceINT) then
    interX=integralOverQYpoint_S_Rec(var,process,incCut,CutParam,QMin_in,Q3,yMin_in,yMax_in,X1,X2,X3,valueMax)&
    +integralOverQYpoint_S_Rec(var,process,incCut,CutParam,Q3,Qmax_in,yMin_in,yMax_in,X3,X4,X5,valueMax)
  else
    interX=valueACB
  end if
end function integralOverQYpoint_S_Rec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!THE MAIN INTERFACE TO CROSS-SECTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! interface for array,s,array,array,array,logical,optional, optional
subroutine MainInterface_AsAAAloo(X,process,s,qT,Q,y,includeCuts,CutParameters)
  integer,intent(in),dimension(1:4)::process                    !the number of process
  real(dp),intent(in)::s                    !Mandelshtam s
  real(dp),intent(in),dimension(1:2)::qT            !(qtMin,qtMax)
  real(dp),intent(in),dimension(1:2)::Q                !(Qmin,Qmax)
  real(dp),intent(in),dimension(1:2)::y                !(ymin,ymax)
  logical,intent(in)::includeCuts                !include cuts
  real(dp),intent(in),dimension(1:4),optional::CutParameters    !(p1,p2,eta1,eta2)

  real(dp)::X

  real(dp),dimension(1:4)::CutParam

  !!! determine cut parameters
  if(includeCuts) then
    if(present(CutParameters)) then
      CutParam=CutParameters
    else
      write(*,*) ErrorString('called includeCuts=true, while CutParameters are undefined',moduleName)
      write(*,*) ErrorString('Evaluation stop',moduleName)
      stop
    end if
  else
    CutParam=(/0d0,0d0,0d0,0d0/)
  end if

  !!!! evaluation
  CallCounter=CallCounter+1

  X=Xsec_IntegrateALL(process,includeCuts,CutParam,s,y(1),y(2),Q(1),Q(2),qT(1),qT(2))

end subroutine MainInterface_AsAAAloo
  
  
subroutine xSec_DY_List(X,process,s,qT,Q,y,includeCuts,CutParameters)
  integer,intent(in),dimension(:,:)::process            !the number of process
  real(dp),intent(in),dimension(:)::s                !Mandelshtam s
  real(dp),intent(in),dimension(:,:)::qT            !(qtMin,qtMax)
  real(dp),intent(in),dimension(:,:)::Q                !(Qmin,Qmax)
  real(dp),intent(in),dimension(:,:)::y                !(ymin,ymax)
  logical,intent(in),dimension(:)::includeCuts        !include cuts
  real(dp),intent(in),dimension(:,:)::CutParameters            !(p1,p2,eta1,eta2)
  real(dp),dimension(:),intent(out)::X
  integer :: i,length

  length=size(s)

  !!! cheking sizes
  if(size(X)/=length) then
    write(*,*) ErrorString('xSec_DY_List: sizes of xSec and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(process,1)/=length) then
    write(*,*) ErrorString('xSec_DY_List: sizes of process and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(qT,1)/=length) then
    write(*,*) ErrorString('xSec_DY_List: sizes of qT and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(y,1)/=length) then
    write(*,*) ErrorString('xSec_DY_List: sizes of y and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(Q,1)/=length) then
    write(*,*) ErrorString('xSec_DY_List: sizes of Q and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(includeCuts)/=length) then
    write(*,*) ErrorString('xSec_DY_List: sizes of includeCuts and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(CutParameters,1)/=length) then
    write(*,*) ErrorString('xSec_DY_List: sizes of CutParameters and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(process,2)/=4) then
    write(*,*) ErrorString('xSec_DY_List: process list must be (:,1:4).',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(qT,2)/=2) then
    write(*,*) ErrorString('xSec_DY_List: qt list must be (:,1:2).',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(y,2)/=2) then
    write(*,*) ErrorString('xSec_DY_List: y list must be (:,1:2).',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(Q,2)/=2) then
    write(*,*) ErrorString('xSec_DY_List: Q list must be (:,1:2).',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if

  CallCounter=CallCounter+length

  !$OMP PARALLEL DO DEFAULT(SHARED)

    do i=1,length
      X(i)=Xsec_IntegrateALL(process(i,1:4),includeCuts(i),CutParameters(i,1:4),s(i),y(i,1),y(i,2),Q(i,1),Q(i,2),qT(i,1),qT(i,2))
    end do
  !$OMP END PARALLEL DO

end subroutine xSec_DY_List
  
subroutine xSec_DY_List_BINLESS(X,process,s,qT,Q,y,includeCuts,CutParameters)
  integer,intent(in),dimension(:,:)::process            !the number of process
  real(dp),intent(in),dimension(:)::s                !Mandelshtam s
  real(dp),intent(in),dimension(:)::qT            !(qtMin,qtMax)
  real(dp),intent(in),dimension(:)::Q                !(Qmin,Qmax)
  real(dp),intent(in),dimension(:)::y                !(ymin,ymax)
  logical,intent(in),dimension(:)::includeCuts        !include cuts
  real(dp),intent(in),dimension(:,:)::CutParameters            !(p1,p2,eta1,eta2)
  real(dp),dimension(:),intent(out)::X

  real(dp),allocatable,dimension(:,:)::vv
  integer :: i,length

  length=size(s)

  !!! cheking sizes
  if(size(X)/=length) then
    write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of xSec and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(process,1)/=length) then
    write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of process and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(qT)/=length) then
    write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of qT and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(y)/=length) then
    write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of y and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(Q)/=length) then
    write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of Q and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(includeCuts)/=length) then
    write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of includeCuts and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(CutParameters,1)/=length) then
    write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of CutParameters and s lists are not equal.',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if
  if(size(process,2)/=4) then
    write(*,*) ErrorString('xSec_DY_List_BINLESS: process list must be (:,1:4).',moduleName)
    write(*,*) ErrorString('Evaluation stop',moduleName)
    stop
  end if

  CallCounter=CallCounter+length

  allocate(vv(1:length,1:7))

  !$OMP PARALLEL DO DEFAULT(SHARED)

    do i=1,length
      vv(i,1:7)=kinematicArray(qt(i),s(i),Q(i),y(i))
      X(i)=xSec(vv(i,1:7),process(i,1:4),includeCuts(i),CutParameters(i,1:4))

    end do
  !$OMP END PARALLEL DO
  deallocate(vv)
end subroutine xSec_DY_List_BINLESS
  
end module TMDX_DY
