!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.0
!
!    This module contains the object for working with grids of optimalTMDs
!       i.e. it stores a function of (real,real,-5:5,1:hadron)
!       The logarithmic scaling in applied for both directtons
!       In principal any such function can be used, but main application is to store (optimal) TMD_OPE results.
!       The grid is Chebishev in directions, with subgrids
!
!                A.Vladimirov (16.10.2025)
!---------------------------------------------------------------------------------------

module aTMDe_optGrid
use aTMDe_Numerics
use aTMDe_interfaces
use aTMDe_IO
implicit none

private

character(len=7),parameter :: moduleName="optGrid"

!!!!!! The object of optimal-TMD grid --------
type, public :: optGrid
  private
  !!!!indicator that grid is ready to use.
  logical::gridReady=.false.
  !!!!!! this parameters needed for generation of messages
  character(:), allocatable::parentName
  integer::outputLevel
  !!!! xRanges & bRanges are the list of values of subgrids for x and b
  real(dp),allocatable::xRanges(:),bRanges(:)
  !!!! number of Subgrids : numbering 1->numXsubgrids
  integer::numXsubgrids,numBsubgrids
  !!!! xNodes & bNodes are the list of values of nodes in the terms of T=[-1,1]
  real(dp),allocatable::xNodes(:),bNodes(:)
  !!!! xNodes & bNodes are the list of factors (-1)^i/2^%., there %=1 for i=0,NUM, and 0 otherwice
  real(dp),allocatable::xNodeFactors(:),bNodeFactors(:)
  !!!! xGridSize & bGridSize are the number of nodes in the subgrids
  !!!! number of nodes is made same for all subgrids, in order to simplify memory operation (store all nodes in single multi-array)
  !!!! : numbering 0->xGridSize
  integer::xGridSize,bGridSize
  !!!! utmost values of the grids
  real(dp)::xMIN,bMIN,bMAX
  !!!! number of hadrons
  integer::numH
  !!!! include Gluon
  logical::withGluon


  !!! parameter of tolerance
  real(dp)::zero=1.d-16

  !!!! xIntervals & bIntervals are the list of (u_{k+1}-u_k)/2 for subgrids
  !!!! xMeans & bMeanss are the list of (u_{k+1}+u_k)/2 for subgrids
  !!!! these are used to speed-up transformation from and to grids (the list are in the trnasformed variables)
  real(dp),allocatable::xIntervals(:),bIntervals(:),xMeans(:),bMeans(:)

  !!!! the main grid variable
  !!!! (subgrid in X,subgrid in B, grid in X, grid in B, flavor(-5:5), hadron)
  real(dp),allocatable::gridMain(:,:,:,:,:,:)
contains
  procedure,public::MakeGrid =>ChGrid_MakeGrid
  procedure,public::Extract  =>ExtractFromGrid
end type

interface optGrid
    procedure :: constructor
end interface optGrid

contains

!!!!! constructor for the optGrid.
!!!!! the initialization happen automatically by reading the constant file.
!!!!! One should specify the path (path), section (moduleLine) and subsection (gridLine) of the constant file where to read
function constructor(path,moduleLine,gridLine,numH_in,withGluon_in,name,outLevel) result(this)
type(optGrid)::this
character(*),intent(in)::path
character(len=5),intent(in)::moduleLine,gridLine
integer,intent(in)::outLevel,numH_in
logical,intent(in)::withGluon_in
character(*),intent(in)::name

integer::i

this%parentName=name
this%outputLevel=outLevel
this%gridReady=.false.

if(this%outputLevel>2) write(*,*) "Initializing optGrid object for "//trim(this%parentName)

OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    !-------------Parameters of grid
    call MoveTO(51,moduleLine)
    call MoveTO(51,gridLine)
    call MoveTO(51,'*p1  ')
    read(51,*) this%numXsubgrids
    allocate(this%xRanges(0:this%numXsubgrids))
    call MoveTO(51,'*p2  ')
    read(51,*) this%xRanges
    call MoveTO(51,'*p3  ')
    read(51,*) this%xGridSize
    call MoveTO(51,'*p4  ')
    read(51,*) this%numBsubgrids
    allocate(this%bRanges(0:this%numBsubgrids))
    call MoveTO(51,'*p5  ')
    read(51,*) this%bRanges
    call MoveTO(51,'*p6  ')
    read(51,*) this%bGridSize
CLOSE (51, STATUS='KEEP')


! !!!processing input parameters
this%xMIN=this%xRanges(0)
this%bMIN=this%bRanges(0)
this%bMAX=this%bRanges(this%numBsubgrids)
this%numH=numH_in

this%withGluon=withGluon_in

!!!! check that the upper limit is zero
if(abs(this%xRanges(this%numXsubgrids)-1._dp)>this%zero) then
  error stop ErrorString('Upper x-range must be 1. Got '//numToStr(this%xRanges(this%numXsubgrids)),this%parentName,moduleName)
end if

!!!!allocation of lists

allocate(this%xNodes(0:this%xGridSize))
allocate(this%bNodes(0:this%bGridSize))
allocate(this%xNodeFactors(0:this%xGridSize))
allocate(this%bNodeFactors(0:this%bGridSize))

allocate(this%xIntervals(1:this%numXsubgrids),this%xMeans(1:this%numXsubgrids))
allocate(this%bIntervals(1:this%numBsubgrids),this%bMeans(1:this%numBsubgrids))

!allocate(this%gridMain(1:this%numXsubgrids,1:this%numBsubgrids,0:this%xGridSize,0:this%bGridSize,-5:5,1:this%numH))

!!!!filing the working variables

this%xIntervals=(log(this%xRanges(1:this%numXsubgrids))-log(this%xRanges(0:this%numXsubgrids-1)))/2._dp
this%xMeans=(log(this%xRanges(1:this%numXsubgrids))+log(this%xRanges(0:this%numXsubgrids-1)))/2._dp

this%bIntervals=(log(this%bRanges(1:this%numBsubgrids))-log(this%bRanges(0:this%numBsubgrids-1)))/2._dp
this%bMeans=(log(this%bRanges(1:this%numBsubgrids))+log(this%bRanges(0:this%numBsubgrids-1)))/2._dp

!!!! node factor is (-1)^i beta_i
this%xNodeFactors=1._dp

do i=0,this%xGridSize
  this%xNodes(i)=cos(i*pi/this%xGridSize)

  if(i==0 .or. i==this%xGridSize) then
    this%xNodeFactors(i)=this%xNodeFactors(i)/2
  end if
  if(mod(i,2)==1) this%xNodeFactors(i)=-this%xNodeFactors(i)

end do

this%bNodeFactors=1._dp

do i=0,this%bGridSize
  this%bNodes(i)=cos(i*pi/this%bGridSize)

  if(i==0 .or. i==this%bGridSize) then
    this%bNodeFactors(i)=this%bNodeFactors(i)/2
  end if
  if(mod(i,2)==1) this%bNodeFactors(i)=-this%bNodeFactors(i)
end do


if(this%outputLevel>2) then
    write(*,*) 'Grid options:'
    write(*,'(A)',advance="no")  ' |  xSubGrid                 ='
    write(*,'(ES10.3)') this%xRanges
    write(*,'(A)',advance='no')  ' |  bSubGrid                 ='
    write(*,'(ES10.3)') this%bRanges
    write(*,'(A,I6,A,I6,A)') ' |  (GridSizeX,GridSizeB)    =(',this%xGridSize,',',this%bGridSize,')'
    write(*,'(A,I3)')   ' |  hadrons to grid           =',this%numH
end if
if(this%outputLevel>2) write(*,*) "Initialization of optGrid object for "//trim(this%parentName)//" complete."

end function constructor

!!!!! values of x computed from the nodes
!!!!! n=subgrid, k=node
pure function XfromNode(this,n,k)
class(optGrid), intent(in)::this
integer,intent(in)::n,k
real(dp)::XfromNode
XfromNode=exp(this%xIntervals(n)*this%xNodes(k)+this%xMeans(n))
end function XfromNode

!!!!! values of b computed from the nodes
!!!!! n=subgrid, k=node
pure function BfromNode(this,n,k)
class(optGrid), intent(in)::this
integer,intent(in)::n,k
real(dp)::BfromNode
BfromNode=exp(this%bIntervals(n)*this%bNodes(k)+this%bMeans(n))
end function BfromNode

!!! this subroutine creates the grid.
!!! F is the function of CxF_compute(x_local,b_local,h,withGluon)
subroutine ChGrid_MakeGrid(this,F)
class(optGrid), intent(inout)::this
  procedure(optTMD_proc)::F
  real(dp):: x_local,b_local
  integer:: iX,iB,h,jX,jB,n
  real(dp)::time1,time2
  !$ real*8::omp_get_wtime

  this%gridReady=.false.

  !!!!!! reallocation of the grid
  if(allocated(this%gridMain)) deallocate(this%gridMain)
  allocate(this%gridMain(1:this%numXsubgrids,1:this%numBsubgrids,0:this%xGridSize,0:this%bGridSize,-5:5,1:this%numH))

  call cpu_time(time1)
  !$ time1=omp_get_wtime()
  if(this%outputlevel>2) write(*,*) 'arTeMiDe.',this%parentName,'.optGrid starts to compute grid.'

  do h=1,this%numH
   !!! !$OMP PARALLEL DO PRIVATE(iB, jX, jB, n, x_local, b_local,h)
   do iX=1,this%numXsubgrids
    do iB=1,this%numBsubgrids
      do jX=0,this%xGridSize
        do jB=0,this%bGridSize
          x_local=XfromNode(this,iX,jX)
          b_local=BfromNode(this,iB,jB)

          this%gridMain(iX,iB,jX,jB,-5:5,h)=F(x_local,b_local,h)

          do n=-5,5
           if(ISNAN(this%gridMain(iX,iB,jX,jB,n,h))) then
            write(*,'(" Function value at (x,b,f,h) =(",F6.5,", ",F6.2,", ",I2,", ",I2,") is computed to NaN")') &
            x_local,b_local,h,n
            end if

            if(abs(this%gridMain(iX,iB,jX,jB,n,h))>10.d6) then
            write(*,'(" Function value at (x,b,f,h) =(",F6.5,", ",F6.2,", ",I2,", ",I2,") is > 10^6")') &
            x_local,b_local,h,n
            write(*,*) "-->",this%gridMain(iX,iB,jX,jB,-5:5,h)
            end if
          end do

          end do
        end do
      end do
    end do
    !!! !$OMP END PARALLEL DO
    if(this%outputLevel>1 .and. this%numH>1) write(*,'(" ",A,": optGrid for hadron ",I3," is done")') this%parentName,h
   end do

  call cpu_time(time2)
  !$ time2=omp_get_wtime()

  if(this%outputlevel>1) then
    if(this%numH>1) then
      write(*,'(" ",A,": optGrids are built  (",I3," x",I3,")x(",I3," x",I3,")  calc.time=",F6.2,"s. ")')&
        this%parentName, this%numXsubgrids,this%xGridSize,this%numBsubgrids,this%bGridSize, time2-time1
    else
      write(*,'(" ",A,": optGrid is built  (",I3," x",I3,")x(",I3," x",I3,")  calc.time=",F6.2,"s. ")')&
        this%parentName, this%numXsubgrids,this%xGridSize,this%numBsubgrids,this%bGridSize, time2-time1
    end if
  end if
  this%gridReady=.true.
    
end subroutine ChGrid_MakeGrid


!!!!! this function interpolates the the grid to the element t
!!!!! it is the baricentric formula
!!!!! f(x)=sum b(i)f(i))/(t-t(i))/sum b(i)/(t-t(i))
pure function interpolateInX(this,t,grid)
class(optGrid), intent(in)::this
real(dp),intent(in),dimension(0:this%xGridSize,-5:5)::grid
real(dp),intent(in)::t
real(dp),dimension(-5:5)::interpolateInX
real(dp),dimension(0:this%xGridSize)::deltaT
integer::i

deltaT=t-this%xNodes
!!! pass though each term, if it is close to 0, then the value of fucntion is given by this node
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
!!!!! just as the previous formula but it operates over list (size of b-nodes) of grids (speed up by factor n)
!!!!! returns the list of interpolation values
pure function interpolateInX_array(this,t,grid)
class(optGrid), intent(in)::this
real(dp),intent(in),dimension(0:this%xGridSize,0:this%bGridSize,-5:5)::grid
real(dp),intent(in)::t
real(dp),dimension(0:this%bGridSize,-5:5)::interpolateInX_array
real(dp),dimension(0:this%xGridSize)::deltaT
integer::i,j

deltaT=t-this%xNodes
!!! pass though each term, if it is close to 0, then the value of fucntion is given by this node
do i=0,this%xGridSize
  if(abs(deltaT(i))<this%zero) then
    interpolateInX_array(0:this%bGridSize,-5:5)=grid(i,0:this%bGridSize,-5:5)
    return
  end if
end do

deltaT=this%xNodeFactors/deltaT
do i=0,this%bGridSize
do j=-5,5
 interpolateInX_array(i,j)=sum(deltaT(:)*grid(:,i,j))
end do
end do
interpolateInX_array=interpolateInX_array/sum(deltaT)

end function interpolateInX_array

!!!!! this function interpolates the the grid to the element t
!!!!! it is the baricentric formula
!!!!! f(x)=sum b(i)f(i))/(t-t(i))/sum b(i)/(t-t(i))
pure function interpolateInB(this,t,grid)
class(optGrid), intent(in)::this
real(dp),intent(in),dimension(0:this%bGridSize,-5:5)::grid
real(dp),intent(in)::t
real(dp),dimension(-5:5)::interpolateInB
real(dp),dimension(0:this%bGridSize)::deltaT
integer::i

deltaT=t-this%bNodes
!!! pass though each term, if it is close to 0, then the value of fucntion is given by this node
do i=0,this%bGridSize
  if(abs(deltaT(i))<this%zero) then
    interpolateInB(-5:5)=grid(i,-5:5)
    return
  end if
end do

deltaT=this%bNodeFactors/deltaT
interpolateInB=matmul(deltaT,grid)/sum(deltaT)
end function interpolateInB

!!!!! this function interpolates the the grid to the element t
!!!!! just as the previous formula but it operates over list (size of x-nodes) of grids (speed up by factor n)
!!!!! returns the list of interpolation values
pure function interpolateInB_array(this,t,grid)
class(optGrid),intent(in)::this
real(dp),intent(in),dimension(0:this%xGridSize,0:this%bGridSize,-5:5)::grid
real(dp),intent(in)::t
real(dp),dimension(0:this%xGridSize,-5:5)::interpolateInB_array
real(dp),dimension(0:this%bGridSize)::deltaT
integer::i,j

deltaT=t-this%bNodes
!!! pass though each term, if it is close to 0, then the value of fucntion is given by this node
do i=0,this%bGridSize
  if(abs(deltaT(i))<this%zero) then
    interpolateInB_array(0:this%xGridSize,-5:5)=grid(0:this%xGridSize,i,-5:5)
    return
  end if
end do

deltaT=this%bNodeFactors/deltaT
do i=0,this%xGridSize
do j=-5,5
 interpolateInB_array(i,j)=sum(deltaT(:)*grid(i,:,j))
end do
end do
interpolateInB_array=interpolateInB_array/sum(deltaT)

end function interpolateInB_array

!!!! return the value of t=[-1,1] from value of x, and number of subgrid
pure function TfromX(this,x,n)
class(optGrid),intent(in)::this
integer,intent(in)::n
real(dp),intent(in)::x
real(dp)::TfromX

TfromX=(log(x)-this%xMeans(n))/this%xIntervals(n)

end function TfromX

!!!! return the value of t=[-1,1] from value of b, and number of subgrid
pure function TfromB(this,b,n)
class(optGrid),intent(in)::this
integer,intent(in)::n
real(dp),intent(in)::b
real(dp)::TfromB

TfromB=(log(b)-this%bMeans(n))/this%bIntervals(n)
end function TfromB

!!!! Interpolation from the grid.
!!!! for values below xMin, above x=1 terminates
!!!! for values below bMin freeze value
!!!! for values above bMax =0
pure function ExtractFromGrid(this,x,bT,h)
class(optGrid),intent(in)::this
  real(dp),intent(in)::x,bT
  integer,intent(in)::h
  real(dp),dimension(-5:5)::ExtractFromGrid

  integer::i,nX,nB
  real(dp)::tX,tB
  real(dp),dimension(0:this%xGridSize,-5:5)::interGrid
  !real(dp),dimension(0:this%bGridSize,-5:5)::interGrid

  if(.not.this%gridReady) then
    error stop ErrorString('attempt to extract from grid while it is not ready',this%parentName,moduleName)
  end if

  !!! checking exeptions
  if(h==0 .or. h>this%numH) then
    error stop ErrorString('the hadron '//numToStr(h)//' is not found in the grid',this%parentName,moduleName)
  end if

  if(x+this%zero<this%XMin) then
   error stop ErrorString('The TMD with x ='//numToStr(x)//' is called. Current grid size is up to '//&
              numToStr(this%XMin)//'. Enlarge boundaries.',this%parentName,moduleName)
  end if
  if(x>1.d0) then
   error stop ErrorString('The TMD with x >1 ('//numToStr(x)//') is called.',this%parentName,moduleName)
  end if
  if(bT<0d0) then
   error stop ErrorString('The TMD with bT <0 ('//numToStr(bT)//') is called.',this%parentName,moduleName)
  end if

  if(x==1.d0) then
    ExtractFromGrid=0._dp
    return
  end if

  !!!! searching for the subgrid in X
  nX=0
  do i=1,this%numXsubgrids
    if(x<this%xRanges(i)) then
      nX=i
      exit
    end if
  end do

  !!!!! for b>bMAX freeze the value
  if(bT>=this%BMAX) then
   tX=TfromX(this,x,nX)
   ExtractFromGrid=interpolateInX(this,tX,this%gridMain(nX,this%numBsubgrids,0:this%xGridSize,0,-5:5,h))

  else
    !!!! searching for the subgrid in B
    !!!! note that nB=0 implies that 0<b<bMIN
    nB=0
    do i=0,this%numBsubgrids
      if(bT<this%bRanges(i)) then
        nB=i
        exit
      end if
    end do

    !!!! case if 0<b<bMIN, return value at bMIN
    if(nB==0) then
      tX=TfromX(this,x,nX)
      ExtractFromGrid=interpolateInX(this,tX,this%gridMain(nX,1,0:this%xGridSize,this%bGridSize,-5:5,h))

    else

      tX=TfromX(this,x,nX)
      tB=TfromB(this,bT,nB)

      !!! first in b then in x
      interGrid=interpolateInB_array(this,tB,this%gridMain(nX,nB,0:this%xGridSize,0:this%bGridSize,-5:5,h))
      ExtractFromGrid=interpolateInX(this,tX,interGrid)
      !!! first in x then in b
      !interGrid=interpolateInX_array(this,tX,this%gridMain(nX,nB,0:this%xGridSize,0:this%bGridSize,-5:5,h))
      !ExtractFromGrid=interpolateInB(this,tB,interGrid)
    end if
  end if
!   do i=-5,5
!     if(ISNAN(ExtractFromGrid(i))) then
!      write(*,'(" Extracted value (x,b,f,h) =(",F6.5,", ",F6.2,", ",I2,", ",I2,") is computed to NaN")') &
!             x,bT,h,i
!       write(*,*) "==>",ExtractFromGrid
!       write(*,*) " "
!       write(*,*) "-->",interGrid(:,i)
!       write(*,*) " "
!       write(*,*) "T->",tX,tB
!       write(*,*) " "
!       write(*,*) "X->",tX-xNodes
!       write(*,*) " "
!       write(*,*) "B->",tB-bNodes
!       stop
!     end if
!     if(abs(ExtractFromGrid(i))>1000000.d0) then
!      write(*,'(" Extracted value (x,b,f,h) =(",F6.5,", ",F6.2,", ",I2,", ",I2,") is computed >10^6")') &
!             x,bT,h,i
!       write(*,*) "==>",ExtractFromGrid
!       write(*,*) " "
!       write(*,*) "-->",interGrid(:,i)
!       write(*,*) " "
!       write(*,*) "T->",tX,tB
!       write(*,*) " "
!       write(*,*) "X->",tX-xNodes
!       write(*,*) " "
!       write(*,*) "B->",tB-bNodes
!       stop
!     end if
!   end do

end function ExtractFromGrid

end module aTMDe_optGrid
