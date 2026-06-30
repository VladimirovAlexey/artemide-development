!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.05
!
!    This is a part of TMDX_DY complex of modules.
!    This module computes the integral over the bin for cross-section of a Drell-Yan reaction (and related sub-processes)
!    The cross-section is evaluated in TMDX_DY_point
!
!    if you use this module please, quote 1706.01473
!
!    ver 3.05: created from v.3.04 (AV, 29.06.2026)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDX_DY_bin
use aTMDe_Numerics
use aTMDe_math
use aTMDe_IO
use aTMDe_ptSpec
use aTMDe_Integration
use EWinput
use TMDX_DY_point

implicit none
private

!!!!!! 1=accurate but slow
!!!!!! 2=fast but not accurate
#define INTEGRATION_MODE 2

!Current version of module
character (len=11),parameter :: moduleName="TMDX-DY-bin"
character (len=5),parameter :: version="v3.05"
!Last appropriate version of constants-file
integer,parameter::inputver=42

real(dp) :: toleranceINT=0.0001d0
real(dp) :: toleranceGEN=0.0000001d0

integer::outputlevel=2
type(Warning_OBJ)::Warning_Handler

!!! number of sections for PT-integral by default
integer::NumPTdefault=6
!!! Maximum size of Q-bin. Larger bins are desected
real::maxQbinSize=30.
!!! Minimal qT, below this number the value is frozen
real(dp)::qTMin_ABS=0.0001d0

!!!--------- bin specifications
!!! below these values the computation will be done for a point
real(dp)::qTBINSIZE_min=0.01_dp
real(dp)::QBINSIZE_min=0.01_dp
real(dp)::yBINSIZE_min=0.0001_dp  !!!! also used for xF


!!!--------- approximate evaluation of qT-bin integration
!!! Idea is based on the fact, that cross-section curve is very smooth vs. qT,
!!! so few-point Chebyshev approximation gives very good result and can be split into bins and integrated exactly
!!! Min Number of Chebyschev nodes in a range (for ranges bigger than 20 GeV, the number is increased automatically)
logical::doPartitioning_byDefault=.false.
integer::NumChNodes=10
real(dp)::MaxQT_range_toPartite=15._dp  !!!!! if the range is bigger than this, it is cut
!!! the matrix of interpolation, and Chebyschev nodes
real(dp),allocatable,dimension(:,:)::ChInterpolationMatrix
real(dp),allocatable,dimension(:)::ChNodes


logical::started=.false.

public:: TMDX_DY_bin_Initialize,TMDX_DY_bin_IsInitialized,TMDX_DY_bin_ResetCounters
public:: Xsec_PTint_Qint_Yint,Xsec_PTspectrum_Qint_Yint

contains

function TMDX_DY_bin_IsInitialized()
  logical::TMDX_DY_bin_IsInitialized
  TMDX_DY_bin_IsInitialized=started
end function TMDX_DY_bin_IsInitialized

!! Initialization of the package
subroutine TMDX_DY_bin_Initialize(file,prefix)
  character(len=*),intent(in)::file
  character(len=*),intent(in),optional::prefix
  character(len=:),allocatable::path
  logical::initRequired
  integer::i,j,FILEver,messageTrigger

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

  !!!! Part A* is about definition of factorization details
  !!!------ factorization parameters are in TMDX_DY_point.

  !!!------ parameters of numerics
  call MoveTO(51,'*B   ')
  call MoveTO(51,'*p1  ')
  read(51,*) toleranceGEN
  call MoveTO(51,'*p2  ')
  read(51,*) toleranceINT
  call MoveTO(51,'*p3  ')
  read(51,*) NumPTdefault
  call MoveTO(51,'*p4  ')
  read(51,*) maxQbinSize
  call MoveTO(51,'*p5  ')
  read(51,*) qTMin_ABS
  call MoveTO(51,'*p6  ')
  read(51,*) doPartitioning_byDefault
  call MoveTO(51,'*p7  ')
  read(51,*) NumChNodes
  call MoveTO(51,'*p8  ')
  read(51,*) MaxQT_range_toPartite

  !!!------ parameters of bins
  call MoveTO(51,'*BIN ')
  call MoveTO(51,'*p1  ')
  read(51,*) qTBINSIZE_min
  call MoveTO(51,'*p2  ')
  read(51,*) QBINSIZE_min
  call MoveTO(51,'*p3  ')
  read(51,*) yBINSIZE_min

  !!!! Part C* and D* are about specifications of factorization implmentations
  !!!------ factorization parameters are in TMDX_DY_point.

  CLOSE (51, STATUS='KEEP')

  if(.not.EWinput_IsInitialized()) then
    if(outputLevel>1) write(*,*) '.. initializing EWinput (from ',moduleName,')'
    if(present(prefix)) then
      call EWinput_Initialize(file,prefix)
    else
      call EWinput_Initialize(file)
    end if
  end if
  if(.not.TMDX_DY_1pt_IsInitialized()) then
    if(outputLevel>1) write(*,*) '.. initializing TMDX_DY_point (from ',moduleName,')'
    if(present(prefix)) then
      call TMDX_DY_1pt_Initialize(file,prefix)
    else
      call TMDX_DY_1pt_Initialize(file)
    end if
  end if

  Warning_Handler=Warning_OBJ(moduleName=moduleName,messageCounter=0,messageTrigger=messageTrigger)

  !!!! the definition of Chebyshev interpolation matrix, which interpolates the qT-sequences
  allocate(ChInterpolationMatrix(0:NumChNodes,0:NumChNodes),ChNodes(0:NumChNodes))
  do i=0,NumChNodes
  ChNodes(i)=cos(i*pi/NumChNodes)
  do j=0,NumChNodes
    ChInterpolationMatrix(i,j)=Cos(i*j*pi/NumChNodes)*2/NumChNodes
    if(i==0 .or. i==NumChNodes) ChInterpolationMatrix(i,j)=ChInterpolationMatrix(i,j)/2
    if(j==0 .or. j==NumChNodes) ChInterpolationMatrix(i,j)=ChInterpolationMatrix(i,j)/2
  end do
  end do

  started=.true.

#if INTEGRATION_MODE==2
    write(*,*)  color('--------------------------------------------------------',c_red)
    write(*,*)  color('----------------------  WARNING!  ----------------------',c_red)
    write(*,*)  color('--   TMDX_DY is in the approximate integration mode   --',c_red)
    write(*,*)  color('--            Faster, but lower precision.            --',c_red)
    write(*,*)  color('--    Switch to default version by changing flag      --',c_red)
    write(*,*)  color('-- INTEGRATION_MODE in TMDX_DY_bin.f90, and recompile --',c_red)
    write(*,*)  color('--------------------------------------------------------',c_red)
#endif

  write(*,*)  color('----- arTeMiDe.TMD_DY_bin '//trim(version)//': .... initialized',c_green)
end subroutine TMDX_DY_bin_Initialize

!!!!! Reset counters here
subroutine TMDX_DY_bin_ResetCounters()
    call Warning_Handler%Reset()
end subroutine TMDX_DY_bin_ResetCounters

!!! check is the process y-symmetric
pure function IsySymmetric(p2)
  logical::IsySymmetric
  integer,intent(in)::p2
  if(p2==1 .or. p2==3 ) then
    IsySymmetric=.true.
  else
    IsySymmetric=.false.
  end if
end function IsySymmetric


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS CALCULATING CROSS-SECTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! ------------------------------------------------------------------------------------------------------------
!!!!! In this module, I compute AVARAGE over bin cross-section,
!!!!! i.e. I integrate over bin, and then divide by the size of the bin
!!!!! It is important, because then the result of computation over a point is essentially the same as computation over a small-bin (i.e. scaling is preserved)
!!!!!
!!!!! ------------------------------------------------------------------------------------------------------------
!!!!! the integration checks/pass the boundaries as following
!!!!! there are bin-boundaries (y1,y2), (Q1,Q2), (qT1,qT2), and the current DYpoint p
!!!!! 1) Integration of qT. if qT1=qT2 = no integration pass to intY_intQ(p)
!!!!!                       else pass to integral over qt with update of DYpoint. intY_intQ(p.updateQT) [it does not change Q,y]
!!!!! 2) Integration of Q. if Q1=Q2 = no integration pass to intY(p)
!!!!!                       else pass to integral over Q with update of DYpoint. intY(p.updateQ) [it does not change y]
!!!!! 3) Integration of Y. if y1=y2 = no integration pass to xSec(p)
!!!!!                       else pass to integral over y with update of DYpoint. xSec(p.updateY)
!!!!! In this way the parameter p always changes. IT IS NOT PARALLEL_SAFE (!!!!); so keep all integrals in a succesive order
!!!!! Also, this structure automatically handle all posibilities of integrations over (y,Q,qT).

!---------------------------------INTEGRATED over Y---------------------------------------------------------------
!!! function for integration over Y
!!! It is crucially important to do the integral over y the first,
!!! because if the integration is done over xF, its boundary depends on current qT and Q.
!!! Note that DYpoint carries information about actual y; while limits are given in "what is needed"
function Xsec_Yint(p,process,incCut,CutParam,yMIN_in,yMAX_in)
  type(DYpoint),intent(inout)::p
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process
  real(dp), intent(in)::yMIN_in,yMAX_in !!!! it is expected that limits are strongly ordered

  real(dp) :: ymin,ymax
  real(dp) :: Xsec_Yint
  real(dp) :: ymin_Check,ymax_Check

  if(isConvergenceLost_inSubmodule()) then
    Xsec_Yint=1d9
    return
  end if

  !!!! ------- NO INTEGRATION PART --------
  !!!! the the boundaries are equal, there is no need to compute the integral
  !!!! in principle here can appear x1,x2>1, but it is nullified at the level of TMD_F
  if(abs(ymax_in-ymin_in)<yBINSIZE_min) then
    !!!! if there is no integration one should set y to its actual value (to be sure, who known what is inside..)
    if(process(1)==2) then
        !!!!! For the integration over xF
        !!!!! The DY point should be recomputed from xF to y, the process is unchanged because it contains Jacobian
        call updateY_DY(p,yFromXF(ymax_in,p))
    else
        !!!!! no integration
        call updateY_DY(p,ymax_in)
    end if
    Xsec_Yint=xSec_DY_1pt(p,process,incCut,CutParam)
    return
  end if

  !!!! ------- WITH INTEGRATION PART --------
  !!! Compute the limits of integration
  !!! in the case process=2 the integral is over xF
  if(process(1)==2) then
    ymin=yFromXF(ymin_in,p)
    ymax=yFromXF(ymax_in,p)
    !!
    !!process(1)=1 !!!! this is important because the actual integration is over y, and process=2 contains Jacobian.
  else
    ymin=ymin_in
    ymax=ymax_in
  end if

  !!!!! the variables x1,x2 could be larger than 1; so there are limiting values of y=\pm log(tau). (note that tau<1)
  !!!!! Note, that in the case of non-exactX1X2 kinematics the tau is different and boundary is a bit larger than 1
  !!!!! in this case, the cross-section is nullified at the level of TMDF_F.
  ymin_Check=log(p%tau)
  ymax_Check=-log(p%tau)

  !!! checking that the integral is inside of physical limits
  if(ymax<ymin_check .or. ymin>ymax_check) then
    Xsec_Yint=0._dp
    return
  end if
  !!! Restricting the limits the cases
  if(ymax > ymax_check) then
    ymax=yMax_check
  end if
  if(ymin < ymin_check) then
    ymin=ymin_check
  end if

  !!!!! In the case of symmetric processes in the symetric limits
  !!!!! I integrate only over the half of integral and multiply by 2
  if(IsySymmetric(process(1)) .and. (ABS(ymax+ymin)<toleranceGEN)) then
    !!! factor 2 is due to symetry condition
#if INTEGRATION_MODE==1
    !!!! slower but accurate
    Xsec_Yint=2*Integrate_SA(integrandOverY,0._dp,ymax,toleranceINT)
#elif INTEGRATION_MODE==2
   !!!! fast but not that accurate
    Xsec_Yint=2*Integrate_G7(integrandOverY,0._dp,ymax)
#endif

  else!!!!! usual integration
#if INTEGRATION_MODE==1
    !!!! slower but accurate
    Xsec_Yint=Integrate_SA(integrandOverY,ymin,ymax,toleranceINT)
#elif INTEGRATION_MODE==2
    !!!! fast but not that accurate
    Xsec_Yint=Integrate_G7(integrandOverY,ymin,ymax)
#endif

  end if

  !!!!! returns average over the bin
  !!!!! IMPORTANT!!!! if the bin goes beyond the physical region, the integration is still weighted by bin-size without restrictions.
  !!!!! It is done to easily handle "complete-y" integration.
  !!!!! ANOTHER IMPORTANT!!! for the case proc1=2 (dxF cross-section) the average over xF-bin is done by \Delta xF.
  !!!!! Then in the limit \Delta xF -> 0, it reproduces the Jacobian, and returns correct cross-section d\sigma/dx_F
  !!!!! which coincides with no-integration path
  Xsec_Yint=Xsec_Yint/(ymax_in-ymin_in)

contains

function integrandOverY(y)
real(dp),intent(in)::y
real(dp)::integrandOverY
call updateY_DY(p,y)

!!!!! BUG CORRECTED: 2.01.2026; Very many thanks to Valentin Moos!
!!!!! importantly, the process is always dy thus process(1)->1, strictly
!!!!! the reason is that all changes in the Jacbian are translated into the changes of limits of integration
integrandOverY=xSec_DY_1pt(p,(/1,process(2),process(3),process(4)/),incCut,CutParam)
end function integrandOverY

end function Xsec_Yint

!---------------------------------INTEGRATED over Y over Q---------------------------------------------------------------
!!!! No need for check over Y they take place within y-integration for each value of Q(!)
!!!!
!!!! to integrate over Q, I use adaptive Simpson.
!!!! However, before it, I check the size of Q-bin, if Q-bin is large I split it into sub-bins (this way helps to avoid problem of missed MZ-peak)
!!!! If the integration is not accurate (i.e. approximate integration mode)
!!!! I separate out Z-peak (+-3GeV), and integrate it with SA, in order not to lose precision too much
!!!!
!!!! Another note, the integral is over dQ while the xSec is over dQ2;
function Xsec_Qint_Yint(p,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in)
  type(DYpoint),intent(inout)::p
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process
  real(dp),intent(in) :: yMin_in,yMax_in  !!!! these varialbes are just passing to intY
  real(dp),intent(in) :: QMin_in,QMax_in  !!!! limits of integration, they are expected t obe strongly ordered
  real(dp):: Xsec_Qint_Yint
  integer::numSec,i
  real(dp)::dQ

  if(isConvergenceLost_inSubmodule()) then
    Xsec_Qint_Yint=1d9
    return
  end if

  !!!! ------- NO INTEGRATION PART --------
  !!!! the the boundaries are equal, there is no need to compute the integral
  if(abs(QMax_in-QMin_in)<QBINSIZE_min) then
    !!!! if there is no integration one should set Q to its actual value (to be sure, who known what is inside..)
    !!!! It is the cross-section d\sigma/dQ
    call updateQ_DY(p,QMax_in)
    Xsec_Qint_Yint=2*p%Q*Xsec_Yint(p,process,incCut,CutParam,yMin_in,yMax_in)
    return
  end if

  !!!! ------- WITH INTEGRATION PART --------

!!!! slower but accurate
#if INTEGRATION_MODE==1
  !!! check how many maxQbins is inside the integration range (+1)
  numSec=INT((Qmax_in-Qmin_in)/maxQbinSize)+1

  !!! if the bin is smaller than maxQbinSize, integrate as is
  if(numSec==1) then

    Xsec_Qint_Yint=Integrate_SA(integrandOverQ,Qmin_in,Qmax_in,toleranceINT)
  else
    !!! else divide to smaler bins and sum the integrals
    dQ=(Qmax_in-Qmin_in)/numSec !!! size of new bins

    Xsec_Qint_Yint=0._dp
    do i=0,numSec-1
      Xsec_Qint_Yint=Xsec_Qint_Yint + &
          Integrate_SA(integrandOverQ,Qmin_in+i*dQ,Qmin_in+(i+1)*dQ,toleranceINT)
    end do
  end if
!!!! fast but not that accurate
#elif INTEGRATION_MODE==2
!!!! in this case I only check for the Z-boson peak
!!!! if it is in the region, I integrate over it specially
  if(Qmin_in<MZ-3 .and. MZ+3<QMax_in) then
    Xsec_Qint_Yint=Integrate_G7(integrandOverQ,Qmin_in,MZ-3)
    Xsec_Qint_Yint=Xsec_Qint_Yint+Integrate_SA(integrandOverQ,MZ-3,MZ+3,toleranceINT)
    Xsec_Qint_Yint=Xsec_Qint_Yint+Integrate_G7(integrandOverQ,MZ+3,Qmax_in)
  else
    Xsec_Qint_Yint=Integrate_G7(integrandOverQ,Qmin_in,Qmax_in)
  end if
#endif

  !!!! returns average over the bin
  Xsec_Qint_Yint=Xsec_Qint_Yint/(Qmax_in-Qmin_in)

contains

function integrandOverQ(Q)
real(dp),intent(in)::Q
real(dp)::integrandOverQ
!!!!! Note that it contains the Jacobian 2Q dQ=dQ^2
call updateQ_DY(p,Q)
integrandOverQ=2*Q*Xsec_Yint(p,process,incCut,CutParam,yMin_in,yMax_in)
end function integrandOverQ

end function Xsec_Qint_Yint

!---------------------------------INTEGRATED over Y over Q over pT-------------------------------------------------------------
!!! No need for check over Y they take a place within y-integration for each value of pT(!)
!!! the checks of PT and Q happen before
!!!!
!!!! to integrate over qT, I use adaptive Simpson (in the pricise version). The integration starts with Num/2 sections.
!!!! If the integration is not accurate (i.e. approximate integration mode) the integration is made by Num-sections (Num = even)
!!!!
!!!! Another note, the integral is over dqT2 while the xSec is over dqT2;
function Xsec_PTint_Qint_Yint_calc(p,process,incCut,CutParam,qTMin,qTMax,Qmin,Qmax,ymin,ymax,Num)
type(DYpoint),intent(inout)::p
logical,intent(in)::incCut
real(dp),dimension(1:4),intent(in)::CutParam
integer,dimension(1:4),intent(in)::process
real(dp), intent(in) :: qTMin,qTMax
real(dp),intent(in):: Qmin,Qmax,ymin,ymax!!!!! there parameters just pass to further integrations
integer,intent(in) :: Num
real(dp):: Xsec_PTint_Qint_Yint_calc

integer::numSec,i
real(dp)::dqT

if(isConvergenceLost_inSubmodule()) then
  Xsec_PTint_Qint_Yint_calc=1d9
  return
end if

!!!! ------- NO INTEGRATION PART --------
!!!! the the boundaries are equal, there is no need to compute the integral
if(abs(qTMax-qTMin)<qTBINSIZE_min) then
    !!!! if there is no integration one should set qT to its actual value (to be sure, who known what is inside..)
    !!!! It is the cross-section d\sigma/dqT
    call updateQT_DY(p,qTMax)
    Xsec_PTint_Qint_Yint_calc=2*p%qT*Xsec_Qint_Yint(p,process,incCut,CutParam,Qmin,Qmax,ymin,ymax)
    return
end if

!!!! ------- WITH INTEGRATION PART --------
!!!! slower but accurate
#if INTEGRATION_MODE==1
!!! since SA starts from 4 sections, and then is going into two subdivision.
!!! So, I should divide Num by 4 to get the same minimal precision in the adaptive case (+1 to secure case of N=2)
  numSec=Num/4+1

  !!! if the bin is smaller than typical range of change, integrate as is
  if(numSec==1) then
    Xsec_PTint_Qint_Yint_calc=Integrate_SA(integrandOverQT,qTMin,qTMax,toleranceINT)
  else
    !!! else divide to smaler bins and sum the integrals
    dqT=(qTMax-qTmin)/numSec !!! size of new bins

    Xsec_PTint_Qint_Yint_calc=0._dp
    do i=0,numSec-1
      Xsec_PTint_Qint_Yint_calc=Xsec_PTint_Qint_Yint_calc + &
          Integrate_SA(integrandOverQT,qTMin+i*dqT,qTMin+(i+1)*dqT,toleranceINT)
    end do
  end if

!!!! fast but not that accurate
#elif INTEGRATION_MODE==2
Xsec_PTint_Qint_Yint_calc=Integrate_SN(integrandOverQT,qTMin,qTMax,Num)
#endif

!!!! returns average over the bin
Xsec_PTint_Qint_Yint_calc=Xsec_PTint_Qint_Yint_calc/(qTMax-qTMin)

contains

function integrandOverQT(qT)
real(dp),intent(in)::qT
real(dp)::integrandOverQT
call updateQT_DY(p,qT)
integrandOverQT=2*p%qT*Xsec_Qint_Yint(p,process,incCut,CutParam,Qmin,Qmax,ymin,ymax)
end function integrandOverQT

end function Xsec_PTint_Qint_Yint_calc

!!!!
!!!!  This routine is entry point to the module
!!!! It checks conditions of the proper bin, creates and sends it to the computation
function Xsec_PTint_Qint_Yint(process,incCut,CutParam,s_in,qt_min_in,qt_max_in,Q_min_in,Q_max_in,ymin_in,ymax_in,Num)
logical,intent(in)::incCut
real(dp),dimension(1:4),intent(in)::CutParam
integer,dimension(1:4),intent(in)::process
integer,intent(in) :: Num
real(dp),intent(in):: ymin_in,ymax_in,Q_min_in,Q_max_in,qt_min_in,qt_max_in,s_in

real(dp):: Xsec_PTint_Qint_Yint
real(dp):: Q_min,Q_max,qt_min,qt_max,s
type(DYpoint)::p

!!!------------------------- checking Q----------
if(Q_min_in<0.9d0) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Q<0.9.')
  write(*,*) "Qmin =",Q_min_in," (Qmin set to 1.GeV)"
  Q_min=1._dp
else
  Q_min=Q_min_in
end if

!!!!!! Comparison of boundaries is made relative to tolerance,
!!!!!!   because they can be very close to each other indicating non-itegration mode
if(Q_min-Q_max_in>toleranceGEN) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Qmax<Qmin. RESULT 0')
  Xsec_PTint_Qint_Yint=0._dp
  return
end if
Q_max=Q_max_in

!!!------------------------- checking S----------
if(s_in<0.9d0) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with s<0.9.')
  write(*,*) "s =",s_in," (s set to Qmin)"
  s=Q_min
else if(s_in<Q_min) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with s<Qmin')
  write(*,*) "s =",s_in," (s set to Qmin)"
  s=Q_min
else
  s=s_in
end if

!!!------------------------- checking PT----------
if(qT_min_in<0.0d0) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with qT<0.')
  write(*,*) "qTmin =",qT_min_in," (qTmin set to 0.GeV)"
  qT_min=qTMin_ABS
else if(qT_min_in<qTMin_ABS) then
    !!!!! if qT is between zero and qTMIN_ABS, freeze to qTMIN_ABS
    !!!!! this is normal situation -- no Warnings
  qT_min=qTMin_ABS
else
  qT_min=qT_min_in
end if

if(qT_min-qT_max_in>toleranceGEN) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with qTmax<qTmin. RESULT 0')
  Xsec_PTint_Qint_Yint=0._dp
  return
end if
qT_max=qT_max_in

!!!------------------------- checking Y----------

if(ymin_in-ymax_in>toleranceGEN) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Ymax<Ymin. RESULT 0')
  Xsec_PTint_Qint_Yint=0._dp
  return
end if


p=DYpoint(qT=qT_min,s=s,Q=(Q_min+Q_max)/2,y=(ymin_in+ymax_in)/2)

Xsec_PTint_Qint_Yint=Xsec_PTint_Qint_Yint_calc(p,process,incCut,CutParam,qT_Min,qT_Max,Q_min,Q_max,ymin_in,ymax_in,Num)

end function Xsec_PTint_Qint_Yint

!---------------------------------INTEGRATED over Y over Q over pT-------------------------------------------------------------
!!!!! This is analog of Xsec_PTint_Qint_Yint but with PT-bins
!!! integration over PT is made by interpolation of the whole range and integrating exactly the Chebyshev interpolation
!!! the bins are presented as list of minimal, and maximal values. All other parameters are same for all bins
function Xsec_PTspectrum_Qint_Yint(process,incCut,CutParam,s_in,qt_min_in,qt_max_in,Q_min_in,Q_max_in,ymin_in,ymax_in)
logical,intent(in)::incCut
real(dp),dimension(1:4),intent(in)::CutParam
integer,dimension(1:4),intent(in)::process
real(dp),intent(in):: ymin_in,ymax_in,Q_min_in,Q_max_in,s_in
real(dp),dimension(1:),intent(in):: qt_min_in,qt_max_in
real(dp),dimension(1:size(qt_min_in)):: Xsec_PTspectrum_Qint_Yint

integer::nBINS,i,j
type(DYpoint)::p
real(dp):: Q_min,Q_max,qt_min(1:size(qt_min_in)),qt_max(1:size(qt_min_in)),s
real(dp)::qT_bin_MAX,qT_bin_MIN,diffAB,sumAB,qTInter
real(dp),dimension(0:NumChNodes)::fAtNodes,vectorR
!!!!! scale for the rescaling of large integrals
real(dp),parameter::scaleL=2.5_dp

if(isConvergenceLost_inSubmodule()) then
  Xsec_PTspectrum_Qint_Yint=1d9
  return
end if

nBINS=size(qt_min_in)
if(size(qt_max_in)/=nBINS) then
  error stop ErrorString('Sizes of qT-min and qT-max lists do not coincide',moduleName)
end if

!!!------------------------- checking Q----------
if(Q_min_in<0.9d0) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Q<0.9.')
  write(*,*) "Qmin =",Q_min_in," (Qmin set to 1.GeV)"
  Q_min=1._dp
else
  Q_min=Q_min_in
end if

if(Q_min-Q_max_in>toleranceGEN) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Qmax<Qmin. RESULT 0')
  Xsec_PTspectrum_Qint_Yint=0._dp
  return
end if
Q_max=Q_max_in

!!!------------------------- checking S----------
if(s_in<0.9d0) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with s<0.9.')
  write(*,*) "s =",s_in," (s set to Qmin)"
  s=Q_min
else if(s_in<Q_min) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with s<Qmin')
  write(*,*) "s =",s_in," (s set to Qmin)"
  s=Q_min
else
  s=s_in
end if

!!!------------------------- checking PT----------
do i=1,nBINS
  if(qT_min_in(i)<0.0d0) then
    call Warning_Handler%WarningRaise('Attempt to compute xSec with qT<0.')
    write(*,*) "qTmin =",qT_min_in(i)," (qTmin set to 0.GeV)"
    qT_min(i)=qTMin_ABS
  else if(qT_min_in(i)<qTMin_ABS) then
    !!!!! if qT is between zero and qTMIN_ABS, freeze to qTMIN_ABS
    !!!!! this is normal situation -- no Warnings
    qT_min(i)=qTMin_ABS
  else
    qT_min(i)=qT_min_in(i)
  end if
  if(qT_min(i)-qT_max_in(i)>toleranceGEN) then
    call Warning_Handler%WarningRaise('Attempt to compute xSec with qTmax<qTmin. RESULT 0')
    Xsec_PTspectrum_Qint_Yint=0._dp
    return
  end if
  !!!!!! It is also importnat to check that bin must be integrated, because PTspectrum integrates all of them
  if(qT_max_in(i)-qT_min(i)<qTBINSIZE_min) then
    call Warning_Handler%WarningRaise('Attempt to compute zero-size bin in PTspectrum. RESULT 0')
    Xsec_PTspectrum_Qint_Yint=0._dp
    return
  end if

  qT_max(i)=qT_max_in(i)
end do


!!!------------------------- checking Y----------

if(ymin_in-ymax_in>toleranceGEN) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Ymax<Ymin. RESULT 0')
  Xsec_PTspectrum_Qint_Yint=0._dp
  return
end if

!!!!!------------------------ here we start the actual computation
!!!!! I have made 2 implementations (1) linear (2) logarithmic
!!!!! determining the size of total qT-range
qT_bin_MIN=minval(qT_min)
qT_bin_MAX=maxval(qT_max)
!!!! makeing intermidiate variables.
!!!---(1) linear
diffAB=(qT_bin_MAX-qT_bin_MIN)/2
sumAB=(qT_bin_MAX+qT_bin_MIN)/2
!!!---(2) logarithmic
!diffAB=(Log(qT_bin_MAX/scaleL+1)-Log(qT_bin_MIN/scaleL+1))/2
!sumAB=(Log(qT_bin_MAX/scaleL+1)+Log(qT_bin_MIN/scaleL+1))/2

!!!!!! Computation of the function at nodes
!!!! longest part computation of the cross-section
do i=0,NumChNodes
  !!!!! ----(1) linear
  qTInter=diffAB*ChNodes(i)+sumAB
  p=DYpoint(qT=qTInter,s=s,Q=(Q_min+Q_max)/2,y=(ymin_in+ymax_in)/2)
  fAtNodes(i)=2*qTInter*Xsec_Qint_Yint(p,process,incCut,CutParam,Q_min,Q_max,ymin_in,ymax_in)

  !!!!! ----(2) logarithm
  !qTInter=scaleL*(exp(diffAB*ChNodes(i)+sumAB)-1)
  !var=kinematicArray(qTInter,s_in,(Q_min+Q_max)/2,(ymin_in+ymax_in)/2)
  !fAtNodes(i)=2*qTInter*(qTInter+scaleL)*Xsec_Qint_Yint(var,process,incCut,CutParam,Q_min,Q_max,ymin_in,ymax_in)
end do

!!!! multiply by the interpolation matrix
fAtNodes=matmul(ChInterpolationMatrix,fAtNodes)
!!!! compute the integral
do i=1,nBINS
   !!!!! ----(1) linear
   vectorR=ChebyshevT_int_array(NumChNodes,(qT_max(i)-sumAB)/diffAB)-ChebyshevT_int_array(NumChNodes,(qT_min(i)-sumAB)/diffAB)
   !!!!! ----(2) logarithm
   !vectorR=ChebyshevT_int_array(NumChNodes,(log(qT_max(i)/scaleL+1)-sumAB)/diffAB) &
   ! -ChebyshevT_int_array(NumChNodes,(log(qT_min(i)/scaleL+1)-sumAB)/diffAB)

  Xsec_PTspectrum_Qint_Yint(i)=diffAB*dot_product(vectorR,fAtNodes)
  !write(*,*) "---->",vectorR

  !!!! do not forget to divide by bin size to get average
  Xsec_PTspectrum_Qint_Yint(i)=Xsec_PTspectrum_Qint_Yint(i)/(qT_max(i)-qT_min(i))
end do

end function Xsec_PTspectrum_Qint_Yint

end module TMDX_DY_bin
