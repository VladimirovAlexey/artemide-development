!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.03
!
!	Contains the module that realizes the Fourier integration with Ogata quadrature
!
!				A.Vladimirov (14.10.2025)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module aTMDe_Ogata
use aTMDe_Numerics
use aTMDe_interfaces
use aTMDe_IO
implicit none

private

!------------------------------------------Tables of zeros of Bessel------------------------------------------------------------
integer,parameter::Nmax=1000
INCLUDE 'Tables/BesselZero1000.f90'

!!!!! I split the qT over runs qT<qTSegmentationBoundary
!!!!! In each segment I have the ogata quadrature with h=hOGATA*hSegmentationWeight
!!!!! It helps to convergen integrals, since h(optimal) ~ qT
integer,parameter::hSegmentationNumber=7
real(dp),dimension(1:hSegmentationNumber),parameter::hSegmentationWeight=(/0.001d0,0.01d0,0.1d0,1d0,2d0,5d0,10d0/)
real(dp),dimension(1:hSegmentationNumber),parameter::qTSegmentationBoundary=(/0.001d0,0.01d0,0.1d0,1d0,10d0,50d0,200d0/)

!!!!!! The object of Ogata integration --------
type, public :: OgataIntegrator
  !!!------ specified upon construction
  !!!!!! this parameters needed for generation of messages
  character(:), allocatable::parentModuleName
  integer::outputLevel

  !!!------ specified by Prepare tables
  !!!!!! the order of the bessel transform used by this integration routine
  integer::N
  !!!!!! general mass scale
  real(dp)::TMDmass,transformationFactor
  !!!!!! numerical parameters of the evaluation
  real(dp)::hOGATA,tolerance
  !!!weights of ogata quadrature
  real(dp),dimension(1:hSegmentationNumber,1:Nmax)::ww
  !!!nodes of ogata quadrature
  real(dp),dimension(1:hSegmentationNumber,1:Nmax)::bb

contains
  !!!!! this is just transformation 1D->1D
  procedure,public:: Transform => Transform_def
  !!!!! this is transformation of function (-5:5) with extra factor M^{2num}/num!/qT^num.
  !!!!! Basically, the transformation of TMD to kT-space
  procedure,public:: TransformTMD => TransformTMD_def
end type

interface OgataIntegrator
    procedure :: constructor
end interface OgataIntegrator

contains

!!!!! constructor for the OgataIntegrator. It automatically prepare the tables with given h
function constructor(parentModuleName,outputLevel,order_in, tolerance_in,hOGATA_in,TMDmass_in) result(this)
type(OgataIntegrator)::this
character(len=*),intent(in)::parentModuleName
integer,intent(in)::outputLevel,order_in
real(dp),intent(in)::tolerance_in,hOGATA_in,TMDmass_in

integer::i,j
real(dp)::xi,hS!=h*hSegmentationWeight

this%parentModuleName=parentModuleName
this%outputLevel=outputLevel
if(this%outputLevel>2) write(*,*) this%parentModuleName,': preparing Ogata tables'

this%TMDmass=TMDmass_in
this%N=order_in
this%tolerance=tolerance_in
this%hOGATA=hOGATA_in

SELECT CASE(this%N)
CASE(0)
    this%transformationFactor=1._dp
CASE(1)
    this%transformationFactor=this%TMDmass**2
CASE(2)
    this%transformationFactor=this%TMDmass**4/2
CASE(3)
    this%transformationFactor=this%TMDmass**6/6
CASE DEFAULT
    write(*,*) ErrorString('aTMDe_Ogata: Ogata transformation with N>3 is not defined. ',this%parentModuleName),this%N
    error stop
END SELECT

do j=1,hSegmentationNumber
do i=1,Nmax

  hS=this%hOGATA*hSegmentationWeight(j)
  xi=JZero(this%N,i)

  !!! if we too far away in xI*hS, the double exponential grow rapidly.
  !!! and for >6, it generates term 10^{300} and exceed the presision

  if(xi*hS>6.d0) then
      this%bb(j,i)=xi*Tanh(piHalf*Sinh(xi*hS))
      this%ww(j,i)=BESSEL_JN(this%N,this%bb(j,i))/xi/(BESSEL_JN(this%N+1,xi)**2)

  else
      this%bb(j,i)=xi*Tanh(piHalf*Sinh(xi*hS))
      this%ww(j,i)=BESSEL_JN(this%N,this%bb(j,i))/xi/(BESSEL_JN(this%N+1,xi)**2)&
      *(pi*xi*hS*Cosh(xi*hS)/(2d0*Cosh(piHalf*Sinh(xi*hS))**2)+Tanh(piHalf*Sinh(xi*hS)))
  end if

end do
end do

if(this%outputLevel>2) write(*,'(A,I4)') ' | Maximum number of nodes    :',Nmax
if(this%outputLevel>1) write(*,*) this%parentModuleName,': Ogata tables (transformation order ',int4ToStr(this%N),') are prepared'
end function constructor

!!!This is the defining module function
!!! It evaluates the integral
!!!  int_0^infty   b^(n+1) db/2  Jn(b qT) F(b)
!!!
!!! or order of hankel transform is defined by this%N,
!!! the function to integrate is "func_1D"
subroutine Transform_def(this,F,qT,res,ISconvergent)
class(OgataIntegrator), intent(inout)::this
procedure(func_1D)::F
real(dp),intent(in)::qT
real(dp),intent(out)::res
logical,intent(out)::ISconvergent

real(dp)::integral,eps,delta
real(dp)::v1,v2,v3,v4
integer::k,j,Nsegment

integral=0.d0
ISconvergent=.true.

v1=1d0
v2=1d0
v3=1d0
v4=1d0

!!! define segment of qT
do j=1,hSegmentationNumber
  if(qT<qTSegmentationBoundary(j)) exit
end do
if(j>hSegmentationNumber) then
  Nsegment=hSegmentationNumber
else
  Nsegment=j
end if

!!! sum over OGATA nodes
do k=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
  eps=this%ww(Nsegment,k)*(this%bb(Nsegment,k)**(this%N+1))*F(this%bb(Nsegment,k)/qT)

  v4=v3
  v3=v2
  v2=v1
  v1=ABS(eps)

  delta=(v1+v2+v3+v4)
  integral=integral+eps

  !!! here we check that residual term is smaller than already collected integral
  !!! also checking the zerothness of the integral. If already collected integral is null it is null
  !!! Here is potential bug. If the first 10 points give zero (whereas some later points do not), the integral will be zero
  if((delta<this%tolerance*abs(integral) .or. abs(integral)<1d-32) .and. k>=10) exit
end do

!!!!! if we run out of maximum number of nodes.
if(k>=Nmax) then
  if(this%outputlevel>0) WRITE(*,*) WarningString('OGATA quadrature diverge. W decaing too slow? ',this%parentModuleName)
  if(this%outputlevel>1) then
    write(*,*) 'Information over the last call ----------'
    write(*,*) 'bt/qT= ',this%bb(Nsegment,Nmax)/qT, 'qT=',qT, '| segmentation zone=',Nsegment,&
        ' ogata h=',this%hOGATA*hSegmentationWeight(Nsegment)
    write(*,*) 'W=',F(this%bb(Nsegment,Nmax)/qT), 'eps/integral =', eps/integral
    write(*,*) 'residual term=',delta, '>',this%tolerance
    write(*,*) '------------------------------------------'
  end if
  ISconvergent=.false.
end if

!!! result is scaled by qT [because the argument of Bessel was scaled bqT-> B]
res=integral/(qT**(this%N+2))

end subroutine Transform_def

!!!This is the defining module function
!!! It evaluates the integral
!!!  int_0^infty   b db/2pi  J_num(b qT) F1  (b/qT)^num M^{2num}/num!
!!!
!!! or order of hankel transform is defined by this%N,
!!! the function to integrate is "func_1D_array5" i.e. (-5:5) function
function TransformTMD_def(this,F,qT)
class(OgataIntegrator), intent(inout)::this
procedure(func_1D_array5)::F
real(dp),intent(in)::qT
real(dp),dimension(-5:5)::TransformTMD_def

real(dp)::integral(-5:5),eps(-5:5)
real(dp)::v1(-5:5),v2(-5:5),v3(-5:5),v4(-5:5),delta(-5:5)
logical:: partDone(-5:5)
integer::k,j,Nsegment

integral=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)

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

do k=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
    eps=this%ww(Nsegment,k)*(this%bb(Nsegment,k)**(this%N+1))*F(this%bb(Nsegment,k)/qT)

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
        if((delta(j)<this%tolerance*ABS(integral(j)) .or. ABS(integral(j))<1d-32) .and. k>=10) partDone(j)=.true.
    end do
    if(partDone(-5).and.partDone(-4).and.partDone(-3).and.partDone(-2).and.partDone(-1)&
        .and.partDone(0).and.partDone(1).and.partDone(2).and.partDone(3).and.partDone(4).and.partDone(5)) exit

end do

if(k>=Nmax) then
    if(this%outputlevel>0) write(*,*) WarningString('OGATA quadrature diverge. TMD decaing too slow? ',this%parentModuleName)
        if(this%outputlevel>2) then
        write(*,*) 'Information over the last call ----------'
        write(*,*) partDone
        write(*,*) 'bt/qT= ',this%bb(Nsegment,Nmax)/qT, 'qT=',qT, '| segmentation zone=',Nsegment,&
            ' ogata h=',this%hOGATA*hSegmentationWeight(Nsegment)
        write(*,*) 'W=',F(this%bb(Nsegment,k)/qT), 'eps/integral =', eps/integral
        write(*,*) 'v1+v2+v3+v4=',v1+v2+v3+v4, '>',this%tolerance*(ABS(integral(1))+ABS(integral(2)))
        write(*,*) '------------------------------------------'
        end if
end if
!!! result is scaled by qT [because the argument of Bessel was scaled bqT-> B]
!!! the extra factor 1/k^n appears due to the defenition of Fourier-trnaform for TMD-n
TransformTMD_def=integral/(qT**(2*this%N+2))*this%transformationFactor

end function TransformTMD_def

end module aTMDe_Ogata
