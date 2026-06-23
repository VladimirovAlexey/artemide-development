!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.04
!
!    This module contains the object for working with kt-grids (depending on Q)
!       i.e. it stores a function of (real,real,real,-5:5,1:hadron)
!       The logarithmic scaling is applied in all directions
!       In principal any such function can be used, but main application is to store TMD_inKT results
!       The grid is Chebyshev along x and kT, and it is Lagrange  along Q
!
!       It stores f=x*TMD, but returns f/x/kT^2, because it is expected that input is multiplied by kT^2 (from Fourier-Levin)
!
!                A.Vladimirov (23.06.2026)
!---------------------------------------------------------------------------------------

module aTMDe_ktGrid
use aTMDe_Numerics
use aTMDe_interfaces
use aTMDe_IO
implicit none

private

character(len=6),parameter :: moduleName="ktGrid"


!!!!!! The object of optimal-TMD grid --------
type, public :: ktGrid
  private
    !!!!indicator that grid is ready to use.
    logical,public::gridReady=.false.
    !!!!!! this parameters needed for generation of messages
    character(:), allocatable::parentName
    integer::outputLevel

    !!!!---------- Variables about the grid in k-space
    !!!! kRanges are the list of values of subgrids for kT
    !!!! Nodes are the list of values of nodes in the terms of T=[-1,1]
    !!!! NodeFactors are the list of factors (-1)^i/2^%., there %=1 for i=0,NUM, and 0 otherwise
    real(dp),allocatable::kRanges(:),kNodes(:),kNodeFactors(:)
    !!!! number of Subgrids
    integer::numKsubgrids
    !!!! number of nodes, it is the same for all subgrids of given class
    integer::kGridSize

    !!!!---------- Variables about the grid in x-space
    !!!! xRanges are the list of values of subgrids for x
    !!!! Nodes are the list of values of nodes in the terms of T=[-1,1]
    !!!! NodeFactors are the list of factors (-1)^i/2^%., there %=1 for i=0,NUM, and 0 otherwise
    real(dp),allocatable::xRanges(:),xNodes(:),xNodeFactors(:)
    !!!! number of Subgrids
    integer::numXsubgrids
    !!!! number of nodes, it is the same for all subgrids of given class
    integer::xGridSize

    !!!!---------- Variables about the grid in Q-space
    !!!! ordinary cubic interpolation over logarithm scale
    !!!! number of nodes in the grid
    integer::QGridSize
    !!!! step size in the logarithmic space
    real(dp)::Qstep
    !!!! nodes in Q, to save the time of computation
    real(dp),allocatable::QNodes(:)

    !!!!---------- General parameters
    !!!! number of hadrons
    integer::numH
    !!!! include Gluon
    logical::withGluon
    !!!! utmost values of the grids
    real(dp)::QMIN,QMAX,kMIN,kMAX,xMIN,lnQMIN,lnQMAX

    real(dp)::zero=10.d-12

    !!!! Intervals are the lists of (u_{k+1}-u_k)/2 for subgrids
    !!!! Means are the lists of (u_{k+1}+u_k)/2 for subgrids
    !!!! these are used to speed-up transformation from and to grids (the lists are in the transformed variables)
    real(dp),allocatable::xIntervals(:),xMeans(:),kIntervals(:),kMeans(:)

    !!!!--------- MAIN GRID VARIABLE
    !!!! (subX,nodeX,subK,nodeK,nodeQ,f,h)
    real(dp),allocatable::mainGRID(:,:,:,:,:,:,:)

contains
  procedure,public::MakeGrid =>ktGrid_MakeGrid
  procedure,public::Extract  =>ExtractFromGrid
end type

interface ktGrid
    procedure :: constructor
end interface ktGrid

contains

!!!!! constructor for the ktGrid.
!!!!! the initialization happens automatically by reading the constant file.
!!!!! One should specify the path (path), section (moduleLine) and subsection (gridLine) of the constant file where to read
function constructor(path,moduleLine,gridLine,numH_in,withGluon_in,name,outLevel) result(this)
type(ktGrid)::this
character(*),intent(in)::path
character(len=5),intent(in)::moduleLine,gridLine
integer,intent(in)::outLevel,numH_in
logical,intent(in)::withGluon_in
character(*),intent(in)::name

integer::i

this%parentName=name
this%outputLevel=outLevel
this%gridReady=.false.

if(this%outputLevel>2) write(*,*) "Initializing ktGrid object for "//trim(this%parentName)

!!!! read input about x and kT-spaces
OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    !-------------Parameters of grid in X
    call MoveTO(51,moduleLine)
    call MoveTO(51,gridLine)
    call MoveTO(51,'*p3  ')
    read(51,*) this%numXsubgrids
    allocate(this%xRanges(0:this%numXsubgrids))
    call MoveTO(51,'*p4  ')
    read(51,*) this%xRanges
    call MoveTO(51,'*p5  ')
    read(51,*) this%xGridSize

    !-------------Parameters of grid in KT
    call MoveTO(51,'*p6  ')
    read(51,*) this%numKsubgrids
    allocate(this%kRanges(0:this%numKsubgrids))
    call MoveTO(51,'*p7  ')
    read(51,*) this%kRanges
    call MoveTO(51,'*p8  ')
    read(51,*) this%kGridSize

    !-------------Parameters of grid in Q
    call MoveTO(51,'*p9  ')
    read(51,*) this%QMIN
    call MoveTO(51,'*p10 ')
    read(51,*) this%QMAX
    call MoveTO(51,'*p11 ')
    read(51,*) this%QGridSize
CLOSE (51, STATUS='KEEP')

this%xMIN=this%xRanges(0)
this%kMIN=this%kRanges(0)
this%kMAX=this%kRanges(this%numKsubgrids)
this%lnQMIN=log(this%QMIN)
this%lnQMAX=log(this%QMAX)

if(abs(this%xRanges(this%numXsubgrids)-1._dp)>this%zero) then
    error stop ErrorString('Upper x-range must be 1. Got '//numToStr(this%xRanges(this%numXsubgrids)),this%parentName,moduleName)
end if

this%numH=numH_in

this%withGluon=withGluon_in

!!!!allocation of lists

allocate(this%xNodes(0:this%xGridSize))
allocate(this%kNodes(0:this%kGridSize))
allocate(this%QNodes(0:this%QGridSize))
allocate(this%xNodeFactors(0:this%xGridSize))
allocate(this%kNodeFactors(0:this%kGridSize))

allocate(this%xIntervals(1:this%numXsubgrids),this%xMeans(1:this%numXsubgrids))
allocate(this%kIntervals(1:this%numKsubgrids),this%kMeans(1:this%numKsubgrids))

!!!!filling the working variables

this%xIntervals=(log(this%xRanges(1:this%numXsubgrids))-log(this%xRanges(0:this%numXsubgrids-1)))/2._dp
this%kIntervals=(log(this%kRanges(1:this%numKsubgrids))-log(this%kRanges(0:this%numKsubgrids-1)))/2._dp
this%xMeans=(log(this%xRanges(1:this%numXsubgrids))+log(this%xRanges(0:this%numXsubgrids-1)))/2._dp
this%kMeans=(log(this%kRanges(1:this%numKsubgrids))+log(this%kRanges(0:this%numKsubgrids-1)))/2._dp

this%xNodeFactors=1._dp

do i=0,this%xGridSize
  this%xNodes(i)=cos(i*pi/this%xGridSize)

  if(i==0 .or. i==this%xGridSize) then
    this%xNodeFactors(i)=this%xNodeFactors(i)/2
  end if
  if(mod(i,2)==1) this%xNodeFactors(i)=-this%xNodeFactors(i)

end do

this%kNodeFactors=1._dp

do i=0,this%kGridSize
  this%kNodes(i)=cos(i*pi/this%kGridSize)

  if(i==0 .or. i==this%kGridSize) then
    this%kNodeFactors(i)=this%kNodeFactors(i)/2
  end if
  if(mod(i,2)==1) this%kNodeFactors(i)=-this%kNodeFactors(i)
end do


this%Qstep=(this%lnQMAX-this%lnQMIN)/this%QGridSize
!!!!! store the values of Q
do i=0,this%QGridSize
    this%QNodes(i)=exp(this%lnQMIN+i*this%Qstep)
end do

!!!!! allocating our huge grid.
allocate(this%mainGRID(&
    1:this%numXsubgrids,&
    0:this%xGridSize,&
    1:this%numKsubgrids,&
    0:this%kGridSize,&
    0:this%QGridSize,&
    -5:5,&
    1:this%numH))

end function constructor

!!!!! values of x computed from the nodes
!!!!! n=subgrid, k=node
pure function XfromNode(this,n,k)
type(ktGrid),intent(in)::this
integer,intent(in)::n,k
real(dp)::XfromNode
XfromNode=exp(this%xIntervals(n)*this%xNodes(k)+this%xMeans(n))
end function XfromNode

!!!!! values of K computed from the nodes
!!!!! n=subgrid, k=node
pure function KfromNode(this,n,k)
type(ktGrid),intent(in)::this
integer,intent(in)::n,k
real(dp)::KfromNode
KfromNode=exp(this%KIntervals(n)*this%KNodes(k)+this%KMeans(n))
end function KfromNode

!!!!! compute and store grid for any TMDgrid_inKT-compatible function
!!!!! In most ocasions this is provided by Fourier_Levin_array in Fourier_Levin module
!!!!! the Grid is TMD*x*kT^2 (factor kT^2 is in the Fourier_Levin_array)
subroutine ktGrid_MakeGrid(this,F)
class(ktGrid), intent(inout)::this
procedure(TMDgrid_inKT)::F
real(dp),dimension(1:this%numKsubgrids,0:this%kGridSize,-5:5):: receivedValues
integer::h,iQ,iX,jX,ff,iK,jK
real(dp)::time1,time2
!$ real*8::omp_get_wtime

call cpu_time(time1)
!$ time1=omp_get_wtime()
if(this%outputlevel>1) write(*,*) 'arTeMiDe.',this%parentName,' starts to compute grid in KT.'

this%gridReady=.false.

!!!!!! reallocation of the grid
if(allocated(this%mainGRID)) deallocate(this%mainGRID)
!!!!! allocating our huge grid.
allocate(this%mainGRID(&
    1:this%numXsubgrids,&
    0:this%xGridSize,&
    1:this%numKsubgrids,&
    0:this%kGridSize,&
    0:this%QGridSize,&
    -5:5,&
    1:this%numH))

do h=1,this%numH
do iQ=0,this%QGridSize
do iX=1,this%numXsubgrids
do jX=0,this%xGridSize
    !!!! requesting the values of Fourier
    receivedValues=F(XfromNode(this,iX,jX),this%QNodes(iQ),h,this%numKsubgrids,this%kGridSize)
    !!!! Here the stored function is multiplied by x (but not by kT^2)
    this%mainGRID(iX,jX,1:this%numKsubgrids,0:this%kGridSize,iQ,-5:5,h)=receivedValues*XfromNode(this,iX,jX)
end do
end do
end do
end do


call cpu_time(time2)
!$ time2=omp_get_wtime()

do h=1,this%numH
do iQ=0,this%QGridSize
do iX=1,this%numXsubgrids
do jX=0,this%xGridSize
do iK=1,this%numKsubgrids
do jK=0,this%kGridSize
do ff=-5,5
    if(ISNAN(this%mainGRID(iX,jX,iK,jK,iQ,ff,h))) then
        ERROR STOP ErrorString("Element of the grid-inKT is NaN. Evaluation STOP",this%parentName,moduleName)
    end if
    if(abs(this%mainGRID(iX,jX,iK,jK,iQ,ff,h))>1.d8) then
        write(*,*) ErrorString("Element of the grid-inKT is greater than 10^8. Evaluation STOP",this%parentName,moduleName)
        write(*,*) "At X=",XfromNode(this,iX,jX)," kT=",KfromNode(this,iK,jK)," Q=",this%QNodes(iQ), "h=",h
        write(*,*) "Value=",this%mainGRID(iX,jX,iK,jK,iQ,:,h)
        ERROR STOP
    end if
end do
end do
end do
end do
end do
end do
end do

if(this%outputlevel>1) then
if(this%numH>1) then
    write(*,&
    '(" ",A,": Grids in KT for ",I3," hadrons are built  (",I3," x",I3,")x(",I3," x",I3,")x",I3,"  calc.time=",F6.2,"s. ")')&
    this%parentName, this%numH, this%numXsubgrids,this%xGridSize,this%numKsubgrids,this%kGridSize,this%QGridSize, time2-time1
else
    write(*,'(" ",A,": Grid in KT is built  (",I3," x",I3,")x(",I3," x",I3,")x",I3,"  calc.time=",F6.2,"s. ")')&
    this%parentName, this%numXsubgrids,this%xGridSize,this%numKsubgrids,this%kGridSize,this%QGridSize, time2-time1
end if
end if

this%gridReady=.true.

end subroutine ktGrid_MakeGrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!  INTERPOLATION PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! return the value of t=[-1,1] from value of x, and number of subgrid
pure function TfromX(this,x,n)
type(ktGrid),intent(in)::this
integer,intent(in)::n
real(dp),intent(in)::x
real(dp)::TfromX

TfromX=(log(x)-this%xMeans(n))/this%xIntervals(n)

end function TfromX

!!!! return the value of t=[-1,1] from value of kT, and number of subgrid
pure function TfromK(this,kT,n)
type(ktGrid),intent(in)::this
integer,intent(in)::n
real(dp),intent(in)::kT
real(dp)::TfromK

TfromK=(log(kT)-this%kMeans(n))/this%kIntervals(n)
end function TfromK

!!!!! this function interpolates the grid to the element t
!!!!! it is the barycentric formula
!!!!! f(x)=sum b(i)f(i))/(t-t(i))/sum b(i)/(t-t(i))
function interpolateInX(this,t,grid)
type(ktGrid),intent(in)::this
real(dp),intent(in),dimension(0:this%xGridSize,-5:5)::grid
real(dp),intent(in)::t
real(dp),dimension(-5:5)::interpolateInX
real(dp),dimension(0:this%xGridSize)::deltaT
integer::i

deltaT=t-this%xNodes
!!! pass through each term, if it is close to 0, then the value of fucntion is given by this node
do i=0,this%xGridSize
    if(abs(deltaT(i))<this%zero) then
      interpolateInX(-5:5)=grid(i,-5:5)
    return
  end if
end do

deltaT=this%xNodeFactors/deltaT

interpolateInX=matmul(deltaT,grid)/sum(deltaT)

end function interpolateInX

!!!!! this function interpolates the the grid to the element t
!!!!! just as the previous formula but it operates over list (size of x-nodes) of grids (speed up by factor n)
!!!!! returns the list of interpolation values
function interpolateInK_array(this,t,grid)
type(ktGrid),intent(in)::this
real(dp),intent(in),dimension(0:this%xGridSize,0:this%kGridSize,-5:5)::grid
real(dp),intent(in)::t
real(dp),dimension(0:this%xGridSize,-5:5)::interpolateInK_array
real(dp),dimension(0:this%kGridSize)::deltaT
integer::i,j

deltaT=t-this%kNodes
!!! pass through each term, if it is close to 0, then the value of fucntion is given by this node
do i=0,this%kGridSize
    if(abs(deltaT(i))<this%zero) then
      interpolateInK_array(0:this%xGridSize,-5:5)=grid(0:this%xGridSize,i,-5:5)
    return
  end if
end do

deltaT=this%kNodeFactors/deltaT
do i=0,this%xGridSize
do j=-5,5
 interpolateInK_array(i,j)=sum(deltaT(:)*grid(i,:,j))
end do
end do
interpolateInK_array=interpolateInK_array/sum(deltaT)

end function interpolateInK_array

!!!! Interpolation from the grid.
!!!! for values below xMin, above x=1 terminates
!!!! for values below kMin freeze value
!!!! for values above kMAX freeze at kMAX node; the kT^2 division at the end  provides the 1/kT^2 tail.
function ExtractFromGrid(this,x,kT,Q,h)
class(ktGrid),intent(in)::this
real(dp),intent(in)::x,kT,Q
integer,intent(in)::h
real(dp),dimension(-5:5)::ExtractFromGrid

integer::i,nX,nK,nQ
real(dp)::tX,tK,kTOdiv

real(dp),dimension(0:this%xGridSize,-5:5)::interGrid
real(dp),dimension(1:4,-5:5)::interQ
real(dp),dimension(1:4)::deltaTQ

if(.not.this%gridReady) then
ERROR STOP ErrorString('attempt to extract from grid while it is not ready',this%parentName,moduleName)
end if

!!! checking exeptions
if(h==0 .or. h>this%numH) then
ERROR STOP ErrorString('the hadron '//numToStr(h)//' is not found in the grid',this%parentName,moduleName)
end if

if(x+this%zero<this%XMin) then
ERROR STOP ErrorString('The TMD with x ='//numToStr(x)//' is called. Current grid size is up to '//&
numToStr(this%XMin)//'. Enlarge boundaries.',this%parentName,moduleName)
end if
if(x>1.d0) then
ERROR STOP ErrorString('The TMD with x >1 ('//numToStr(x)//') is called.',this%parentName,moduleName)
end if
if(kT<0d0) then
ERROR STOP ErrorString('The TMD with kT <0 ('//numToStr(kT)//') is called.',this%parentName,moduleName)
end if
if(Q<=this%QMIN) then
ERROR STOP ErrorString('The TMD with Q <=QMIN ('//numToStr(Q)//'<'//numToStr(this%QMIN)//') is called.',&
            this%parentName,moduleName)
end if
if(Q>=this%QMAX) then
ERROR STOP ErrorString('The TMD with Q >=QMAX ('//numToStr(Q)//">"//numToStr(this%QMAX)//') is called.',&
            this%parentName,moduleName)
end if

if(x==1.d0) then
ExtractFromGrid=0._dp
return
end if

!!!! determining Q-node
!! The clamping of nQ to [1, QGridSize-2] is silently ensuring 4-point Lagrange always has two valid neighbors on each side.
nQ=int((log(Q)-this%lnQMIN)/this%Qstep)
if(nQ==0) then
    nQ=1
else if(nQ==this%QGridSize-1) then
    nQ=this%QGridSize-2
else if(nQ>=this%QGridSize) then
    ERROR STOP ErrorString('The TMD called outside of the grid',this%parentName,moduleName)
end if

!!!! searching for the subgrid in X
nX=0
do i=1,this%numXsubgrids
    if(x<this%xRanges(i)) then
        nX=i
        exit
    end if
end do

!!!!! Reminder::  mainGRID(1:numXsubgrids,0:xGridSize,1:numKsubgrids,0:kGridSize,0:QGridSize,-5:5,1:numH)

!!!!! for kT>kTMax extract the frozen value (it will be devided by kT^2 in the end, which partially correct the behavior)
if(kT>=this%kMAX) then
    tX=TfromX(this,x,nX)
    interQ(1,:)=interpolateInX(this,tX,this%mainGRID(nX,0:this%xGridSize,this%numKsubgrids,0,nQ-1,-5:5,h))
    interQ(2,:)=interpolateInX(this,tX,this%mainGRID(nX,0:this%xGridSize,this%numKsubgrids,0,nQ-0,-5:5,h))
    interQ(3,:)=interpolateInX(this,tX,this%mainGRID(nX,0:this%xGridSize,this%numKsubgrids,0,nQ+1,-5:5,h))
    interQ(4,:)=interpolateInX(this,tX,this%mainGRID(nX,0:this%xGridSize,this%numKsubgrids,0,nQ+2,-5:5,h))

    kTOdiv=kT
else
    !!!! searching for the subgrid in kT
    !!!! note that nK=0 implies that 0<kt<kMIN
    nK=0
    do i=0,this%numKsubgrids
        if(kT<this%kRanges(i)) then
        nK=i
        exit
        end if
    end do

    !!!! case if 0<k<kMIN, return value at kMIN
    if(nK==0) then
        tX=TfromX(this,x,nX)

        interQ(1,:)=interpolateInX(this,tX,this%mainGRID(nX,0:this%xGridSize,1,this%kGridSize,nQ-1,-5:5,h))
        interQ(2,:)=interpolateInX(this,tX,this%mainGRID(nX,0:this%xGridSize,1,this%kGridSize,nQ-0,-5:5,h))
        interQ(3,:)=interpolateInX(this,tX,this%mainGRID(nX,0:this%xGridSize,1,this%kGridSize,nQ+1,-5:5,h))
        interQ(4,:)=interpolateInX(this,tX,this%mainGRID(nX,0:this%xGridSize,1,this%kGridSize,nQ+2,-5:5,h))

        kTOdiv=this%kMIN
    else

        tX=TfromX(this,x,nX)
        tK=TfromK(this,kT,nK)

        !!! first in kT then in x
        interGrid=interpolateInK_array(this,tK,this%mainGRID(nX,0:this%xGridSize,nK,0:this%kGridSize,nQ-1,-5:5,h))
        interQ(1,:)=interpolateInX(this,tX,interGrid)

        interGrid=interpolateInK_array(this,tK,this%mainGRID(nX,0:this%xGridSize,nK,0:this%kGridSize,nQ-0,-5:5,h))
        interQ(2,:)=interpolateInX(this,tX,interGrid)

        interGrid=interpolateInK_array(this,tK,this%mainGRID(nX,0:this%xGridSize,nK,0:this%kGridSize,nQ+1,-5:5,h))
        interQ(3,:)=interpolateInX(this,tX,interGrid)

        interGrid=interpolateInK_array(this,tK,this%mainGRID(nX,0:this%xGridSize,nK,0:this%kGridSize,nQ+2,-5:5,h))
        interQ(4,:)=interpolateInX(this,tX,interGrid)

        kTOdiv=kT
    end if
end if

!!!!! barycentric interpolation with Lagrange polynomials

deltaTQ=Q-this%QNodes(nQ-1:nQ+2)

!!! pass through each term, if it is close to 0, then the value of function is given by this node
if(abs(deltaTQ(1))<this%zero) then
      ExtractFromGrid=interQ(1,:)
else if(abs(deltaTQ(2))<this%zero) then
      ExtractFromGrid=interQ(2,:)
else if(abs(deltaTQ(3))<this%zero) then
      ExtractFromGrid=interQ(3,:)
else if(abs(deltaTQ(4))<this%zero) then
      ExtractFromGrid=interQ(4,:)
else
!!!!! interpolation
!!!!! these are barycentric weights of Lagrange interpolation
deltaTQ=(/-1,3,-3,1/)/deltaTQ/6._dp
ExtractFromGrid=matmul(deltaTQ,interQ)/sum(deltaTQ)
end if

ExtractFromGrid=ExtractFromGrid/x/kTOdiv**2

end function ExtractFromGrid

end module aTMDe_ktGrid
