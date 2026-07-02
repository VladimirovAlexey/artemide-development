!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.05
!
!    This is a part of TMDX_SIDIS complex of modules.
!    This module computes the integral over the bin for cross-section of SIDIS reaction (and related sub-processes)
!    The cross-section is evaluated in TMDX_SIDIS_point
!
!    if you use this module please, quote 1912.06532
!
!    ver 3.05: created from v.3.04 (AV, 02.07.2026)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDX_SIDIS_bin
use aTMDe_Numerics
use aTMDe_math
use aTMDe_IO
use aTMDe_ptSpec
use aTMDe_Integration
use TMDX_SIDIS_point

implicit none
private

!!!!!! 1=accurate but slow
!!!!!! 2=fast but not accurate
#define INTEGRATION_MODE 2

!Current version of module
character (len=14),parameter :: moduleName="TMDX-SIDIS-bin"
character (len=5),parameter :: version="v3.05"
!Last appropriate version of constants-file
integer,parameter::inputver=43

real(dp) :: toleranceINT=0.0001d0
real(dp) :: toleranceGEN=0.0000001d0

integer::outputlevel=2
type(Warning_OBJ)::Warning_Handler

logical::started=.false.

!!! number of sections for PT-integral by default
real(dp)::ptMIN_ABS=0.00001d0

!!!--------- bin specifications
!!! below these values the computation will be done for a point
real(dp)::pTBINSIZE_min=0.01_dp
real(dp)::QBINSIZE_min=0.01_dp
real(dp)::xBINSIZE_min=0.0001_dp
real(dp)::zBINSIZE_min=0.0001_dp

!!!--------- approximate evaluation of pT-bin integration
!!! Idea is based on the fact, that cross-section curve is very smooth vs. pT,
!!! so few-point Chebyshev approximation gives very good result and can be split into bins and integrated exactly
logical::doPartitioning_byDefault=.false.
integer::NumChNodes=10
real(dp)::MaxQT_range_toPartite=18._dp  !!!!! if the range is bigger than this, it is cut
!!! the matrix of interpolation, and Chebyschev nodes
real(dp),allocatable,dimension(:,:)::ChInterpolationMatrix
real(dp),allocatable,dimension(:)::ChNodes


public::TMDX_SIDIS_bin_Initialize,TMDX_SIDIS_bin_IsInitialized,TMDX_SIDIS_bin_ResetCounters

public::Xsec_pTspectrum_Qint_Xint_Zint,Xsec_PTint_Qint_Xint_Zint



contains

function TMDX_SIDIS_bin_IsInitialized()
    logical::TMDX_SIDIS_bin_IsInitialized
    TMDX_SIDIS_bin_IsInitialized=started
end function TMDX_SIDIS_bin_IsInitialized

!! Initialization of the package
subroutine TMDX_SIDIS_bin_Initialize(file,prefix)
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

    !!! Check the file version
    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")

    call MoveTO(51,'*0   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) FILEver
    if(FILEver<inputver) then
        CLOSE (51, STATUS='KEEP')
        error stop ErrorString('const-file version is too old. Update with artemide.setup.',moduleName)
    end if

    !!! Fill the message system
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger

    !!! Go to TMDX-SIDIS section
    call MoveTO(51,'*10  ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        CLOSE (51, STATUS='KEEP')
        return
    end if

    !!!------ parameters of numerics
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) toleranceGEN
    call MoveTO(51,'*p2  ')
    read(51,*) toleranceINT
    call MoveTO(51,'*p4  ')
    read(51,*) ptMIN_ABS
    call MoveTO(51,'*p5  ')
    read(51,*) doPartitioning_byDefault
    call MoveTO(51,'*p6  ')
    read(51,*) NumChNodes
    call MoveTO(51,'*p7  ')
    read(51,*) MaxQT_range_toPartite

      !!!------ parameters of bins
    call MoveTO(51,'*BIN ')
    call MoveTO(51,'*p1  ')
    read(51,*) pTBINSIZE_min
    call MoveTO(51,'*p2  ')
    read(51,*) QBINSIZE_min
    call MoveTO(51,'*p3  ')
    read(51,*) xBINSIZE_min
    call MoveTO(51,'*p4  ')
    read(51,*) zBINSIZE_min


    CLOSE (51, STATUS='KEEP')

    if(.not.TMDX_SIDIS_1pt_IsInitialized()) then
    if(outputLevel>1) write(*,*) '.. initializing TMDX_SIDIS_point (from ',moduleName,')'
    if(present(prefix)) then
      call TMDX_SIDIS_1pt_Initialize(file,prefix)
    else
      call TMDX_SIDIS_1pt_Initialize(file)
    end if
   end if

    Warning_Handler=Warning_OBJ(moduleName=moduleName,messageCounter=0,messageTrigger=messageTrigger)

    !!!! the definition of Chebyshev interpolation matrix, which interpolates the pT cross-section
    allocate(ChInterpolationMatrix(0:NumChNodes,0:NumChNodes),ChNodes(0:NumChNodes))
    do i=0,NumChNodes
    ChNodes(i)=cos(i*pi/NumChNodes)
    do j=0,NumChNodes
        ChInterpolationMatrix(i,j)=Cos(i*j*pi/NumChNodes)*2/NumChNodes
        if(i==0 .or. i==NumChNodes) ChInterpolationMatrix(i,j)=ChInterpolationMatrix(i,j)/2
        if(j==0 .or. j==NumChNodes) ChInterpolationMatrix(i,j)=ChInterpolationMatrix(i,j)/2
    end do
    end do

#if INTEGRATION_MODE==2
    write(*,*)  color('--------------------------------------------------------',c_red)
    write(*,*)  color('----------------------  WARNING!  ----------------------',c_red)
    write(*,*)  color('-- TMDX_SIDIS is in the approximate integration mode  --',c_red)
    write(*,*)  color('--            Faster, but lower precision.            --',c_red)
    write(*,*)  color('--    Switch to default version by changing flag      --',c_red)
    write(*,*)  color('-- INTEGRATION_MODE in TMDX_SIDIS_bin.f90, recompile  --',c_red)
    write(*,*)  color('--------------------------------------------------------',c_red)
#endif

    started=.true.
    write(*,*)  color('----- arTeMiDe.TMD_SIDIS '//trim(version)//': .... initialized',c_green)
end subroutine TMDX_SIDIS_bin_Initialize

!!!!! Reset counters here
subroutine TMDX_SIDIS_bin_ResetCounters()
    call Warning_Handler%Reset()
end subroutine TMDX_SIDIS_bin_ResetCounters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CUTS RELATED FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The cuts are ymin<y<ymax, Wmin<W2<Wmax

!!! checks the lower bound of x against cut constraints and returns the most restrictive lower bound
!!! xmin is the input lower bound; returns max(xmin, cut-induced lower bounds)
pure function xMinWithCuts(xmin,p,cutParam)
type(SIDISpoint),intent(in)::p
real(dp),dimension(1:4),intent(in)::cutParam
real(dp),intent(in)::xmin
real(dp)::x1,x2,xMinWithCuts

x1=p%Q2/cutParam(2)/p%sM2
x2=p%Q2/(p%Q2+cutParam(4)-p%M2target)

!x1=p%Q2/cutParam(2)/(p%sM2-(0.938)**2)
!x2=p%Q2/(p%Q2+cutParam(4)-(0.938)**2)

xMinWithCuts=max(xmin,x1,x2)
end function xMinWithCuts
  
!!! checks the upper bound of x against cut constraints and returns the most restrictive upper bound
!!! xmax is the input upper bound; returns min(xmax, cut-induced upper bounds)
pure function xMaxWithCuts(xmax,p,cutParam)
type(SIDISpoint),intent(in)::p
real(dp),dimension(1:4),intent(in)::cutParam
real(dp),intent(in)::xmax
real(dp)::x1,x2,xMaxWithCuts
x1=p%Q2/cutParam(1)/p%sM2
x2=p%Q2/(p%Q2+cutParam(3)-p%M2target)

!x1=p%Q2/cutParam(1)/(p%sM2-(0.938)**2)
!x2=p%Q2/(p%Q2+cutParam(3)-(0.938)**2)

xMaxWithCuts=min(xmax,x1,x2)
end function xMaxWithCuts
  
!!! checks the lower bound of Q against cut constraints and returns the most restrictive lower bound
!!! Qmin is the input lower bound; xmin is used to evaluate cut-induced lower bounds on Q
pure function QMinWithCuts(xmin,Qmin,p,cutParam)
type(SIDISpoint),intent(in)::p
real(dp),dimension(1:4),intent(in)::cutParam
real(dp),intent(in)::xmin,Qmin
real(dp)::Q1,Q2,QMinWithCuts

Q1=sqrt(xmin*cutParam(1)*p%sM2)
Q2=sqrt(xmin*(cutParam(3)-p%M2target)/(1d0-xmin))

!Q1=sqrt(xmin*cutParam(1)*(p%sM2-(0.938)**2))
!Q2=sqrt(xmin*(cutParam(3)-(0.938)**2)/(1d0-xmin))

QMinWithCuts=max(Qmin,Q1,Q2)
end function QMinWithCuts
  
!!! checks the upper bound of Q against cut constraints and returns the most restrictive upper bound
!!! Qmax is the input upper bound; xmax is used to evaluate cut-induced upper bounds on Q
pure function QMaxWithCuts(xmax,Qmax,p,cutParam)
type(SIDISpoint),intent(in)::p
real(dp),dimension(1:4),intent(in)::cutParam
real(dp),intent(in)::xmax,Qmax
real(dp)::Q1,Q2,QMaxWithCuts

Q1=Sqrt(xmax*cutParam(2)*p%sM2)
Q2=sqrt(xmax*(cutParam(4)-p%M2target)/(1d0-xmax))

!Q1=Sqrt(xmax*cutParam(2)*(p%sM2-(0.938)**2))
!Q2=sqrt(xmax*(cutParam(4)-(0.938)**2)/(1d0-xmax))

QMaxWithCuts=min(Qmax,Q1,Q2)
end function QMaxWithCuts
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS CALCULATING CROSS-SECTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! ------------------------------------------------------------------------------------------------------------
!!!!! In this module, I compute the AVERAGE of the cross-section over the bin,
!!!!! i.e. I integrate over the bin, and then divide by the size of the bin.
!!!!! It is important, because then the result of a computation at a point is essentially the same as a computation over a small bin (smooth limit is preserved).
!!!!!
!!!!! ------------------------------------------------------------------------------------------------------------
!!!!! the integration checks/passes the boundaries as follows
!!!!! there are bin-boundaries (z1,z2), (x1,x2), (Q1,Q2), (pT1,pT2), and the current SIDISpoint p
!!!!! 1) Integration of pPerp. if pT1=pT2, no integration, pass to intZ_intX_intQ(p)
!!!!!                       else pass to integral over pPerp with update of SIDISpoint. intZ_intX_intQ(p.updatePperp) [it does not change Q,x,z]
!!!!! 2) Integration of Q. if Q1=Q2, no integration, pass to intZ_intX(p)
!!!!!                       else pass to integral over Q with update of SIDISpoint. intZ_intX(p.updateQ) [it does not change x,z]
!!!!! 3) Integration of X. if x1=x2, no integration, pass to intZ(p)
!!!!!                       else pass to integral over x with update of SIDISpoint. intZ(p.updateX) [it does not change z]
!!!!! 4) Integration of Z. if z1=z2, no integration, pass to xSec(p)
!!!!!                       else pass to integral over z with update of SIDISpoint. xSec(p.updateZ)
!!!!! In this way the parameter p always changes. IT IS NOT PARALLEL_SAFE (!!!!); so keep all integrals in successive order.
!!!!! Also, this structure automatically handles all possibilities of integrations over (z,x,Q,pT).

    
!---------------------------------INTEGRATED over Z--------------------------------------------------------------
!!! function for integration over Z
function Xsec_Zint(p,process,zMIN_in,zMAX_in)
type(SIDISpoint)::p
real(dp),intent(in) :: zMin_in,zMax_in
integer,dimension(1:4),intent(in)::process

real(dp) :: Xsec_Zint
real(dp)::zMin,zMax

if(isConvergenceLost_inSubmodule()) then
  Xsec_Zint=1d9
  return
end if


!!!!! check the limits
if(zmax_in>1.d0) then
  zmax=1._dp
else
  zmax=zmax_in
end if
zmin=zmin_in

!!!! ------- NO INTEGRATION PART --------
!!!! the boundaries are equal, there is no need to compute the integral
if(abs(zmax-zmin)<zBINSIZE_min) then
  !!!!! no integration
  call updateZ_SIDIS(p,zmax)
  Xsec_Zint=xSec_SIDIS_1pt(p,process)
  return
end if

!!!! ------- WITH INTEGRATION PART --------
#if INTEGRATION_MODE==1
  !!!! slower but accurate
  Xsec_Zint=Integrate_SA(integrandOverZ,zMin,zMax,toleranceINT)
#elif INTEGRATION_MODE==2
  !!!! fast but not that accurate
  Xsec_Zint=Integrate_G7(integrandOverZ,zMin,zMax)
#endif

!!!!! returns average over the bin
Xsec_Zint=Xsec_Zint/(zmax_in-zmin_in)

contains

function integrandOverZ(z)
real(dp),intent(in)::z
real(dp)::integrandOverZ
call updateZ_SIDIS(p,z)
integrandOverZ=xSec_SIDIS_1pt(p,process)
end function integrandOverZ

end function Xsec_Zint
  
!---------------------------------INTEGRATED over X---------------------------------------------------------------

!!! function for integration over X
!!! If proc1=3, the SIDIS point actually contains y instead of x.
!!! it is a curious fact that until this point the factor gamma2 stored in SIDISpoint is incorrect (however it is corrected within TMDX_SIDIS_1pt).
function Xsec_Xint_Zint(p,process,zMIN_in,zMAX_in,Xmin_in,Xmax_in,doCut,Cuts)
type(SIDISpoint)::p
logical,intent(in)::doCut
real(dp),dimension(1:4),intent(in)::Cuts
real(dp),intent(in):: xmin_in,xmax_in,zMIN_in,zMAX_in
integer,dimension(1:4)::process
real(dp) :: Xsec_Xint_Zint

real(dp)::xmin, xmax

if(isConvergenceLost_inSubmodule()) then
  Xsec_Xint_Zint=1d9
  return
end if

!!! ----- check the input variables
!!! in the case process=3, the point contains y, which is to be transformed to x
!!! evaluate the corresponding x bounds from y
if(process(1)==3) then
  xmin=XfromY(xmin_in,p)
  xmax=XfromY(xmax_in,p)
else
  xmin=xmin_in
  xmax=xmax_in
end if

!! in case of cut we determine the recut values
if(doCut) then
  xmin=xMinWithCuts(xmin,p,Cuts)
  xmax=xMaxWithCuts(xmax,p,Cuts)
end if

!! check the limits and fix it to 1
if(xmax>1.d0) xmax=1._dp

!!!! ------- NO INTEGRATION PART --------
!!!! the boundaries are equal, there is no need to compute the integral
if(abs(xmax-xmin)<xBINSIZE_min) then
  !!!!! no integration
  call updateX_SIDIS(p,xmax)
  Xsec_Xint_Zint=Xsec_Zint(p,process,zMIN_in,zMAX_in)
  return
end if

!!!! It may happen that computation of cuts, invert limit, in this case it means that the integral is 0
if(xmin>xmax) then
  Xsec_Xint_Zint=0._dp
  return
end if

!!!! ------- WITH INTEGRATION PART --------
#if INTEGRATION_MODE==1
    !!!! slower but accurate
    Xsec_Xint_Zint=Integrate_SA(integrandOverX,xMin,xMax,toleranceINT)
#elif INTEGRATION_MODE==2
    !!!! fast but not that accurate
    Xsec_Xint_Zint=Integrate_G7(integrandOverX,xMin,xMax)
#endif


!!!!! returns average over the bin
!!!!! IMPORTANT!!!! if the bin goes beyond the physical region, the integration is still weighted by bin-size without restrictions.
!!!!! ANOTHER IMPORTANT!!! for the case proc1=3 (dy cross-section) the average over y-bin is done by \Delta y.
!!!!! Then in the limit \Delta y -> 0, it reproduces the Jacobian, and returns correct cross-section d\sigma/dy
!!!!! which coincides with no-integration path
Xsec_Xint_Zint=Xsec_Xint_Zint/(xmax_in-xmin_in)

contains

function integrandOverX(x)
real(dp),intent(in)::x
real(dp)::integrandOverX
call updateX_SIDIS(p,x)
!!!!! process(1) is forced to 1: actual integration is over x; all Jacobian effects are absorbed into the limits computed from XfromY above.
integrandOverX=Xsec_Zint(p,(/1,process(2),process(3),process(4)/),zMin_in,zMax_in)
end function integrandOverX

end function Xsec_Xint_Zint
  

!---------------------------------INTEGRATED over Q--------------------------------------------------------------

!!! function for integration over Q
!!! If proc1=2, the SIDIS point actually contains y instead of Q.
!!! it is a curious fact that until this point the factor gamma2 stored in SIDISpoint is incorrect (however it is corrected within TMDX_SIDIS_1pt).
!!!! Another note, the integral is over dQ while the xSec is over dQ2;
function Xsec_Qint_Xint_Zint(p,process,zMIN_in,zMAX_in,Xmin_in,Xmax_in,Qmin_in,Qmax_in,doCut,Cuts)
type(SIDISpoint)::p
logical,intent(in)::doCut
real(dp),dimension(1:4),intent(in)::Cuts
real(dp),intent(in):: xmin_in,xmax_in,zMIN_in,zMAX_in,Qmin_in,Qmax_in
integer,dimension(1:4)::process
real(dp) :: Xsec_Qint_Xint_Zint

real(dp)::Qmin, Qmax

if(isConvergenceLost_inSubmodule()) then
  Xsec_Qint_Xint_Zint=1d9
  return
end if

!!! check input parameters
!!! in the case process=2, evaluate the corresponding Q bounds from y
if(process(1)==2) then
  Qmin=QfromY(Qmin_in,p)
  Qmax=QfromY(Qmax_in,p)
else
  Qmin=Qmin_in
  Qmax=Qmax_in
end if

!! in case of cut we determine the recut values
if(doCut) then
  Qmin=QMinWithCuts(Xmin_in,Qmin,p,Cuts)
  Qmax=QMaxWithCuts(Xmax_in,Qmax,p,Cuts)
end if

!!!! ------- NO INTEGRATION PART --------
!!!! the boundaries are equal, there is no need to compute the integral
if(abs(Qmax-Qmin)<QBINSIZE_min) then
  !!!!! no integration
  call updateQ_SIDIS(p,Qmax)
  !!!! for proc1=1: smooth limit of <integral 2Q dQ / DeltaQ> = 2Q * sigma_1
  !!!! for proc1=2: smooth limit of <integral 2Q sigma_1 dQ / DeltaY> = sigma_2 (no 2Q factor)
  if(process(1)==2) then
    Xsec_Qint_Xint_Zint=Xsec_Xint_Zint(p,process,zMIN_in,zMAX_in,Xmin_in,Xmax_in,doCut,Cuts)
  else
    Xsec_Qint_Xint_Zint=2*p%Q*Xsec_Xint_Zint(p,process,zMIN_in,zMAX_in,Xmin_in,Xmax_in,doCut,Cuts)
  end if
  return
end if

!!!! It may happen that computation of cuts, invert limit, in this case it means that the integral is 0
if(Qmin>Qmax) then
  Xsec_Qint_Xint_Zint=0._dp
  return
end if

!!!! ------- WITH INTEGRATION PART --------

!!!! proc1 is forced to 1 only when proc1=2; for proc1=3 the y->x transformation is already handled in Xsec_Xint_Zint
#if INTEGRATION_MODE==1
    !!!! slower but accurate
    if(process(1)==2) then
      Xsec_Qint_Xint_Zint=Integrate_SA(integrandOverQ_proc2,QMin,QMax,toleranceINT)
    else
      Xsec_Qint_Xint_Zint=Integrate_SA(integrandOverQ,QMin,QMax,toleranceINT)
    end if
#elif INTEGRATION_MODE==2
    !!!! fast but not that accurate
    if(process(1)==2) then
      Xsec_Qint_Xint_Zint=Integrate_G7(integrandOverQ_proc2,QMin,QMax)
    else
      Xsec_Qint_Xint_Zint=Integrate_G7(integrandOverQ,QMin,QMax)
    end if
#endif


!!!!! returns average over the bin
!!!!! IMPORTANT!!!! if the bin goes beyond the physical region, the integration is still weighted by bin-size without restrictions.
!!!!! ANOTHER IMPORTANT!!! for the case proc1=2 (dy cross-section) the average over y-bin is done by \Delta y.
!!!!! Then in the limit \Delta y -> 0, it reproduces the Jacobian, and returns correct cross-section d\sigma/dy
!!!!! which coincides with the no-integration path
Xsec_Qint_Xint_Zint=Xsec_Qint_Xint_Zint/(Qmax_in-Qmin_in)

contains

function integrandOverQ(Q)
real(dp),intent(in)::Q
real(dp)::integrandOverQ
call updateQ_SIDIS(p,Q)
integrandOverQ=2*Q*Xsec_Xint_Zint(p,process,zMin_in,zMax_in,Xmin_in,Xmax_in,doCut,Cuts)
end function integrandOverQ

function integrandOverQ_proc2(Q)
real(dp),intent(in)::Q
real(dp)::integrandOverQ_proc2
call updateQ_SIDIS(p,Q)
!!!!! process(1) is forced to 1: actual integration is over Q (converted from y); all Jacobian effects (2Q factor, limits) are handled here.
integrandOverQ_proc2=2*Q*Xsec_Xint_Zint(p,(/1,process(2),process(3),process(4)/),zMin_in,zMax_in,Xmin_in,Xmax_in,doCut,Cuts)
end function integrandOverQ_proc2

end function Xsec_Qint_Xint_Zint

!---------------------------------INTEGRATED over pT over Q over X over Z (calc)-------------------------------------------------------------
!!!! to integrate over pT, I use adaptive Simpson (in the precise version). The integration starts with Num/4+1 sections.
!!!! If the integration is not accurate (i.e. approximate integration mode) the integration is made by Num sections (Num = even)
!!!!
!!!! Another note, the integrand is over dpT (= dpPerp) while the xSec is differential in dpT^2 (hence the factor 2*pT);
function Xsec_PTint_Qint_Xint_Zint_calc(p,process,zMin,zMax,xMin,xMax,Qmin,Qmax,ptMin,ptMax,doCut,Cuts,Num)
type(SIDISpoint)::p
real(dp),dimension(1:4),intent(in) :: Cuts
logical,intent(in)::doCut
integer,dimension(1:4),intent(in)::process
real(dp),intent(in) :: Qmin,Qmax,xMin,xMax,zMin,zMax,ptMax,ptMin
integer,intent(in)::Num
real(dp) :: Xsec_PTint_Qint_Xint_Zint_calc
integer::numSec,i
real(dp)::dqT

if(isConvergenceLost_inSubmodule()) then
  Xsec_PTint_Qint_Xint_Zint_calc=1d9
  return
end if

!!!! ------- NO INTEGRATION PART --------
!!!! the boundaries are equal, there is no need to compute the integral
if(abs(pTMax-pTMin)<pTBINSIZE_min) then
    !!!! if there is no integration one should set pPerp to its actual value (to be sure what state p is in)
    !!!! It is the cross-section d\sigma/dpT^2, multiplied by 2*pPerp (the Jacobian dpT^2 = 2*pT*dpT)
    call updatePperp_SIDIS(p,pTMax)
    Xsec_PTint_Qint_Xint_Zint_calc=2*p%pPerp*Xsec_Qint_Xint_Zint(p,process,zMIN,zMAX,Xmin,Xmax,Qmin,Qmax,doCut,Cuts)
    return
end if

!!!! ------- WITH INTEGRATION PART --------
!!!! slower but accurate
#if INTEGRATION_MODE==1
!!! SA starts from numSec initial sections and then subdivides by bisection.
!!! So, Num/4+1 initial sections give roughly the same minimal precision as Num sections in SN (+1 to handle N=2)
  numSec=Num/4+1

  !!! if the bin is smaller than typical range of change, integrate as is
  if(numSec==1) then
    Xsec_PTint_Qint_Xint_Zint_calc=Integrate_SA(integrandOverPT,pTMin,pTMax,toleranceINT)
  else
    !!! else divide to smaler bins and sum the integrals
    dqT=(pTMax-pTmin)/numSec !!! size of new bins

    Xsec_PTint_Qint_Xint_Zint_calc=0._dp
    do i=0,numSec-1
      Xsec_PTint_Qint_Xint_Zint_calc=Xsec_PTint_Qint_Xint_Zint_calc + &
          Integrate_SA(integrandOverPT,pTMin+i*dqT,pTMin+(i+1)*dqT,toleranceINT)
    end do
  end if

!!!! fast but not that accurate
#elif INTEGRATION_MODE==2
Xsec_PTint_Qint_Xint_Zint_calc=Integrate_SN(integrandOverPT,pTMin,pTMax,Num)
#endif

!!!! returns average over the bin
Xsec_PTint_Qint_Xint_Zint_calc=Xsec_PTint_Qint_Xint_Zint_calc/(pTMax-pTMin)

contains

function integrandOverPT(pT)
real(dp),intent(in)::pT
real(dp)::integrandOverPT
call updatePperp_SIDIS(p,pT)
integrandOverPT=2*p%pPerp*Xsec_Qint_Xint_Zint(p,process,zMIN,zMAX,Xmin,Xmax,Qmin,Qmax,doCut,Cuts)
end function integrandOverPT

end function Xsec_PTint_Qint_Xint_Zint_calc


!!!!
!!!!  This routine is entry point to the module
!!!! It checks conditions of the proper bin, creates and sends it to the computation
function Xsec_PTint_Qint_Xint_Zint(process,s_in,zMin,zMax,xMin,xMax,Qmin_in,Qmax_in,ptMin_in,ptMax_in,doCut,Cuts,Num,m1,m2)
real(dp),dimension(1:4),intent(in) :: Cuts
logical,intent(in)::doCut
integer,dimension(1:4),intent(in)::process
real(dp),intent(in) :: Qmin_in,Qmax_in,xMin,xMax,zMin,zMax,ptMax_in,ptMin_in,s_in,m1,m2
integer,intent(in)::Num

real(dp) :: Xsec_PTint_Qint_Xint_Zint
real(dp)::Qmin,Qmax,ptMin,ptMax,s
type(SIDISpoint)::p

!!!------------------------- checking Q----------
if(Qmin_in<0.9d0) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Q<0.9.')
  write(*,*) "Qmin =",Qmin_in," (Qmin set to 1.GeV)"
  Qmin=1._dp
else
  Qmin=Qmin_in
end if

!!!!!! Comparison of boundaries is made relative to tolerance,
!!!!!!   because they can be very close to each other indicating non-integration mode
if(Qmin-Qmax_in>toleranceGEN) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Qmax<Qmin. RESULT 0')
  Xsec_PTint_Qint_Xint_Zint=0._dp
  return
end if
Qmax=Qmax_in

!!!------------------------- checking S----------
if(s_in<0.9d0) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with s<0.9.')
  write(*,*) "s =",s_in," (s set to Qmin)"
  s=Qmin
else if(s_in<Qmin) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with s<Qmin')
  write(*,*) "s =",s_in," (s set to Qmin)"
  s=Qmin
else
  s=s_in
end if

!!!------------------------- checking PT----------
if(ptMin_in<0.0d0) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with pPerp<0.')
  write(*,*) "pPerp_min =",ptMin_in," (pPerp_min set to 0.GeV)"
  ptMin=ptMIN_ABS
else if(ptMin_in<ptMIN_ABS) then
    !!!!! if pT is between zero and ptMIN_ABS, freeze to ptMIN_ABS
    !!!!! this is normal situation -- no Warnings
  ptMin=ptMIN_ABS
else
  ptMin=pTmin_in
end if

if(pTmin-pTmax_in>toleranceGEN) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with pTmax<pTmin. RESULT 0')
  Xsec_PTint_Qint_Xint_Zint=0._dp
  return
end if
pTmax=pTmax_in

!!!------------------------- checking Z----------
if(zmin-zmax>toleranceGEN) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Zmax<Zmin. RESULT 0')
  Xsec_PTint_Qint_Xint_Zint=0._dp
  return
end if

!!!------------------------- checking X----------
if(xmin-xmax>toleranceGEN) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Xmax<Xmin. RESULT 0')
  Xsec_PTint_Qint_Xint_Zint=0._dp
  return
end if


p=SIDISpoint(pPerp=ptMin,s=s,Q=(Qmin+Qmax)/2,x=(xMax+xMin)/2,z=(zMax+zMin)/2,Mtarget=m1,Mproduct=m2)


Xsec_PTint_Qint_Xint_Zint=Xsec_PTint_Qint_Xint_Zint_calc(p,process,zMin,zMax,xMin,xMax,Qmin,Qmax,ptMin,ptMax,doCut,Cuts,Num)

end function Xsec_PTint_Qint_Xint_Zint

!---------------------------------INTEGRATED over pT (and Z and Q and X) -------------------------------------------------------------
!!!!! This is analog of Xsec_PTint_Qint_Xint_Zint but evaluates multiple pT-bins simultaneously
!!! integration over pT is made by Chebyshev interpolation of the whole range, then integrating the interpolant exactly
!!! the bins are provided as lists of minimum and maximum values. All other parameters are the same for all bins
function Xsec_pTspectrum_Qint_Xint_Zint(process,s_in,zMin,zMax,xMin,xMax,Qmin_in,Qmax_in,ptMin_in,ptMax_in,doCut,Cuts,m1,m2)
logical,intent(in)::doCut
real(dp),dimension(1:4),intent(in) :: Cuts
integer,dimension(1:4),intent(in):: process
real(dp),intent(in) :: s_in,Qmin_in,Qmax_in,xMin,xMax,zMin,zMax, m1, m2
real(dp), dimension(1:),intent(in) :: ptMin_in, ptMax_in

real(dp), dimension(1:size(ptMin_in)) :: Xsec_pTspectrum_Qint_Xint_Zint
type(SIDISpoint)::p
real(dp) :: s,Qmin,Qmax

integer::nBINS,i
real(dp):: pt_min(1:size(ptMin_in)),pt_max(1:size(ptMin_in))
real(dp)::pT_bin_MAX,pT_bin_MIN,diffAB,sumAB,pTInter
real(dp),dimension(0:NumChNodes)::fAtNodes,vectorR


if(isConvergenceLost_inSubmodule()) then
  Xsec_pTspectrum_Qint_Xint_Zint=1d9
  return
end if

nBINS=size(ptMin_in)
if(size(ptMax_in)/=nBINS) then
  error stop ErrorString('Sizes of pT-min and pT-max lists do not coincide',moduleName)
end if

!!!------------------------- checking pT  ------------------
do i=1,nBINS
  if(ptMin_in(i)<0.0d0) then
    call Warning_Handler%WarningRaise('Attempt to compute xSec with pT<0.')
    write(*,*) "pTmin =",ptMin_in(i)," (pTmin set to 0.GeV)"
    pt_min(i)=ptMIN_ABS
  else if(ptMin_in(i)<ptMIN_ABS) then
    !!!!! if pT is between zero and ptMIN_ABS, freeze to ptMIN_ABS
    !!!!! this is normal situation -- no Warnings
    pt_min(i)=ptMIN_ABS
  else
    pt_min(i)=ptMin_in(i)
  end if
  if(ptMax_in(i)<pt_min(i)) then
    call Warning_Handler%WarningRaise('Attempt to compute xSec with pTmax<pTmin. RESULT 0')
    Xsec_pTspectrum_Qint_Xint_Zint=0._dp
    return
  end if
  if(ptMax_in(i)-pt_min(i)<pTBINSIZE_min) then
      call Warning_Handler%WarningRaise('Attempt to compute zero-size bin in PTspectrum. RESULT 0')
      Xsec_pTspectrum_Qint_Xint_Zint=0._dp
      return
  end if
  pt_max(i)=ptMax_in(i)
end do

!!!------------------------- checking Q----------
if(Qmin_in<0.9d0) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Q<0.9.')
  write(*,*) "Qmin =",Qmin_in," (Qmin set to 1.GeV)"
  Qmin=1._dp
else
  Qmin=Qmin_in
end if

!!!!!! Comparison of boundaries is made relative to tolerance,
!!!!!!   because they can be very close to each other indicating non-integration mode
if(Qmin-Qmax_in>toleranceGEN) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Qmax<Qmin. RESULT 0')
  Xsec_pTspectrum_Qint_Xint_Zint=0._dp
  return
end if
Qmax=Qmax_in

!!!------------------------- checking S----------
if(s_in<0.9d0) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with s<0.9.')
  write(*,*) "s =",s_in," (s set to Qmin)"
  s=Qmin
else if(s_in<Qmin) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with s<Qmin')
  write(*,*) "s =",s_in," (s set to Qmin)"
  s=Qmin
else
  s=s_in
end if

!!!------------------------- checking Z----------
if(zmin-zmax>toleranceGEN) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Zmax<Zmin. RESULT 0')
  Xsec_pTspectrum_Qint_Xint_Zint=0._dp
  return
end if

!!!------------------------- checking X----------
if(xmin-xmax>toleranceGEN) then
  call Warning_Handler%WarningRaise('Attempt to compute xSec with Xmax<Xmin. RESULT 0')
  Xsec_pTspectrum_Qint_Xint_Zint=0._dp
  return
end if


!!!!!------------------------ here we start the actual computation
!!!! determining the total pT range
pT_bin_MIN=minval(pT_min)
pT_bin_MAX=maxval(pT_max)

! write(*,*) pT_bin_MAX
!
!!!! making intermediate variables.
!!!---(1) linear
diffAB=(pT_bin_MAX-pT_bin_MIN)/2
sumAB=(pT_bin_MAX+pT_bin_MIN)/2

!!!!!! Computation of the function at nodes
!!!! longest part: computation of the cross-section at Chebyshev nodes

! write(*,*) "Number of Chebyshev nodes", NumChNodes


do i=0,NumChNodes
  !!!!! ----(1) linear
  pTInter=diffAB*ChNodes(i)+sumAB
  p=SIDISpoint(pPerp=pTInter,s=s,Q=(Qmin+Qmax)/2,x=(xMax+xMin)/2,z=(zMax+zMin)/2,Mtarget=m1,Mproduct=m2)

  fAtNodes(i)=2*pTInter*Xsec_Qint_Xint_Zint(p,process,zMIN,zMAX,Xmin,Xmax,Qmin,Qmax,doCut,Cuts)
end do

!!!! multiply by the interpolation matrix
fAtNodes=matmul(ChInterpolationMatrix,fAtNodes)
!!!! compute the integral
do i=1,nBINS
   !!!!! ----(1) linear
   vectorR=ChebyshevT_int_array(NumChNodes,(pT_max(i)-sumAB)/diffAB)-ChebyshevT_int_array(NumChNodes,(pT_min(i)-sumAB)/diffAB)

   Xsec_pTspectrum_Qint_Xint_Zint(i)=diffAB*dot_product(vectorR,fAtNodes)
   !write(*,*) "---->",vectorR

   !!!! do not forget to divide by bin size to get average
  Xsec_pTspectrum_Qint_Xint_Zint(i)=Xsec_pTspectrum_Qint_Xint_Zint(i)/(pT_max(i)-pT_min(i))
end do

end function Xsec_pTspectrum_Qint_Xint_Zint

end module TMDX_SIDIS_bin
