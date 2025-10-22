!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This file contains the definition of 1D Lagrange grid, which is filled upon request !!!
!!!                             8.10.2025                                               !!!
!!!                                 A.Vladimirov                                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module aTMDe_xGrid
use aTMDe_Numerics
use aTMDe_interfaces
use aTMDe_IO
use EWinput
implicit none

private

!!!!!!!! parameter to check the absence of the value
real(dp),parameter::zero=0.00000000001_dp

type, public::Xgrid !!!!! name of abstract class
    private
     !!!!!! this parameters needed for generation of messages
    character(:), allocatable::parentName
    integer::outputLevel
    type(Warning_OBJ)::Warning_Handler
    !!!!!! this is function to grid
    procedure(strFUNC), pointer, nopass :: FtoGrid
    !!!!!! the process number
    integer,dimension(1:3)::process0
    !!!!!! parameters of the grid
    !!!!!! size of grids
    integer::N1,N2,N3,N4
    !!!!!! minValues and maxValues of grids
    real(dp)::v1MIN,v2MIN,v3MIN,v4MIN
    real(dp)::v1MAX,v2MAX,v3MAX,v4MAX
    !!!!!! extra parameters pf grids
    real(dp)::p11,p21,p31,p41

    real(dp),allocatable::mainGrid(:,:,:,:)

    !!!! counters
    integer::pointsSet=0,pointsCalled=0
contains
    procedure,public::getValue =>getValue_this
    procedure,public::reset =>reset_this
    procedure,private::storePoint
    procedure,private::getValueAtNode
end type

interface Xgrid
    procedure :: constructor
end interface Xgrid

contains

!!!!!!!!!----------------------------------------------------------------------
!!!!!     transformations for arguments to the space of grid [0,N]
!!!!! the trnaformation must satify v=vMin <=>t=0; v=vMax <=>t=N
!!!!! For the moment I use the same transformation for all directions. Namely logarithmic grid
!!!!! t= N* log(v/vMin)/log(vMax/vMin)
!!!!! So, i defined par1=vMin and par2= N/log(vMax/vMin)
pure function VtoT(v,par1,par2)
real(dp),intent(in)::v,par1,par2
real(dp)::VtoT
VtoT=log(v/par1)*par2
end function VtoT

pure function TtoV(t,par1,par2)
real(dp),intent(in)::t,par1,par2
real(dp)::TtoV
TtoV=par1*exp(t/par2)
end function TtoV

!!!!! transformation for qT direction (it is different because qT can be 0) It is quadratice
!!!!! t= N*sqrt((t-vMin)/(vMin-vMax))
!!!!! So, i defined par1=vMin and par2= N/sqrt(vMax-vMin)
pure function VtoT2(v,par1,par2)
real(dp),intent(in)::v,par1,par2
real(dp)::VtoT2
VtoT2=par2*sqrt(v-par1)
end function VtoT2

pure function TtoV2(t,par1,par2)
real(dp),intent(in)::t,par1,par2
real(dp)::TtoV2
TtoV2=(t/par2)**2+par1
end function TtoV2

!!!!!!!!!----------------------------------------------------------------------
!!!!!! Weight function
!!!!!! this function specifies the weight depending on the process
!!!!!! such that grid saves F*W, and this should be a smooth function
function weigthFunction(Q2,qT,x1,x2,process0)
real(dp)::weigthFunction
real(dp),intent(in)::Q2,qT,x1,x2
integer,dimension(1:3),intent(in)::process0

SELECT CASE(process0(3))
CASE(2,3)!!!!!!!! Z-boson peak
    weigthFunction=x1*x2*(qT+1)**2*((Q2-MZ2)**2+GammaZ2*MZ2)/(Q2*Q2)
CASE(4,5,6)!!!!!!!! W-boson peak
    weigthFunction=x1*x2*(qT+1)**2*((Q2-MW2)**2+GammaW2*MW2)/(Q2*Q2)
CASE DEFAULT
    weigthFunction=x1*x2*(qT+1)**2
end select

end function weigthFunction

!!!!! constructor for the Xgrid.
!!!!! one should specify process0, min,max,N for each variable
function constructor(F,proc_in,v1min,v1max,v1N,v2min,v2max,v2N,v3min,v3max,v3N,v4min,v4max,v4N,name,outLevel) result(this)
type(Xgrid)::this
procedure(strFUNC)::F
integer,dimension(1:3),intent(in)::proc_in
integer,intent(in)::v1N,v2N,v3N,v4N
real(dp),intent(in)::v1Min,v1Max,v2Min,v2Max,v3Min,v3Max,v4Min,v4Max
character(*),intent(in)::name
integer::outLevel

this%parentName=name
this%outputLevel=outLevel
if(this%outputLevel>1) write(*,*) trim(this%parentName)//": Initializing the structure-function grid for process ("// &
numToStr(proc_in(1))//","//int4ToStr(proc_in(2))//","//int4ToStr(proc_in(3))//")"

this%Warning_Handler=Warning_OBJ(moduleName=this%parentName,messageCounter=0,messageTrigger=3)

this%pointsSet=0
this%pointsCalled=0

this%FtoGrid=>F
this%process0=proc_in

this%v1MIN=v1MIN
this%v1MAX=v1MAX
this%N1=v1N

this%v2MIN=v2MIN
this%v2MAX=v2MAX
this%N2=v2N

this%v3MIN=v3MIN
this%v3MAX=v3MAX
this%N3=v3N

this%v4MIN=v4MIN
this%v4MAX=v4MAX
this%N4=v4N

this%p11=v1N/Log(v1Max/v1Min)
!this%p21=v2N/Log(v2Max/v2Min)
this%p21=v2N/sqrt(v2Max-v2Min)
this%p31=v3N/Log(v3Max/v3Min)
this%p41=v4N/Log(v4Max/v4Min)

allocate(this%mainGrid(0:this%N1,0:this%N2,0:this%N3,0:this%N3),source=0._dp)

end function constructor

!!!!!!! reset the grid such that it will be computed for anew
subroutine reset_this(this)
class(Xgrid), intent(inout)::this

call this%Warning_Handler%reset()

if(allocated(this%mainGrid)) deallocate(this%mainGrid)
allocate(this%mainGrid(0:this%N1,0:this%N2,0:this%N3,0:this%N3),source=0._dp)

if(this%outputLevel>2) write(*,*) trim(this%parentName)//": grid reset for process ("// &
numToStr(this%process0(1))//","//int4ToStr(this%process0(2))//","//int4ToStr(this%process0(3))//")"
if(this%outputLevel>2) write(*,*) trim(this%parentName)//": nodes used   : "//numToStr(this%pointsSet)
if(this%outputLevel>2) write(*,*) trim(this%parentName)//": nodes called : "//numToStr(this%pointsCalled)

this%pointsSet=0
this%pointsCalled=0

end subroutine reset_this

!!!!!!! computes the value of function at node, n and store it in the MainGrid
!!!!!!! the stored function is multipliedby weigthFunction(x)
subroutine storePoint(this,n1,n2,n3,n4)
class(Xgrid), intent(inout)::this
integer,intent(in)::n1,n2,n3,n4

real*8::v1,v2,v3,v4

this%pointsSet=this%pointsSet+1

v1=TtoV(real(n1,dp),this%v1Min,this%p11)
v2=TtoV2(real(n2,dp),this%v2Min,this%p21)
v3=TtoV(real(n3,dp),this%v3Min,this%p31)
v4=TtoV(real(n4,dp),this%v4Min,this%p41)

this%mainGrid(n1,n2,n3,n4)=this%FtoGrid(v1,v2,v3,v4,this%process0)*weigthFunction(v1,v2,v3,v4,this%process0)
end subroutine storePoint

!!!!!!! Intents to extract the value of function at the node. If it is not setup compute it and stores
function getValueAtNode(this,n1,n2,n3,n4)
class(Xgrid), intent(inout)::this
real(dp)::getValueAtNode
integer,intent(in)::n1,n2,n3,n4

this%pointsCalled=this%pointsCalled+1
!write(*,*) "enter:",n1,n2,n3,n4
!!!! check that function is stored comparing to the default value
if(abs(this%MainGrid(n1,n2,n3,n4))<zero) then
    !write(*,*) "store:",n1,n2,n3,n4
    call this%storePoint(n1,n2,n3,n4)
end if
!!!! return value from the grid.
getValueAtNode=this%MainGrid(n1,n2,n3,n4)

end function getValueAtNode

!!!!!!! Computes the interpolation from the grid
!!!!!!! using the Neuvile algorithm for qubic Legandre interpolation
!!!!!!! it is assumed that distance between node is =1
function getValue_this(this,v1,v2,v3,v4)
class(Xgrid), intent(inout)::this
real(dp)::getValue_this
real(dp),intent(in)::v1,v2,v3,v4

real(dp)::t1,t2,t3,t4
integer::n1,n2,n3,n4,i1,i2,i3,i4
real(dp)::subGrid4D(0:3,0:3,0:3,0:3),subGrid3D(0:3,0:3,0:3),subGrid2D(0:3,0:3),subGrid1D(0:3)

!!!!! check tha point inside of the grid
if(v1<this%v1MIN .or. v2<this%v2MIN .or. v3<this%v3MIN .or. v4<this%v4MIN &
    .or. v1>this%v1MAX .or. v2>this%v2MAX .or. v3>this%v3MAX .or. v4>this%v4MAX) then
    if(v1<this%v1MIN .or. v1>this%v1MAX) call this%Warning_Handler%WarningRaise("Point Q2 outside of specified grid !")
    if(v2<this%v2MIN .or. v2>this%v2MAX) call this%Warning_Handler%WarningRaise("Point qT outside of specified grid !")
    if(v3<this%v3MIN .or. v3>this%v3MAX) call this%Warning_Handler%WarningRaise("Point x1 outside of specified grid !")
    if(v4<this%v4MIN .or. v4>this%v4MAX) call this%Warning_Handler%WarningRaise("Point x2 outside of specified grid !")


    if(this%outputLevel>2) then
        write(*,*) "  Requested point:", v1,v2,v3,v4
        write(*,*) "  Min value      :",this%v1MIN,this%v2MIN,this%v3MIN,this%v4MIN
        write(*,*) "  Max value      :",this%v1MAX,this%v2MAX,this%v3MAX,this%v4MAX
    end if
    getValue_this=this%FtoGrid(v1,v2,v3,v4,this%process0)
end if

!write(*,*) "ask for:",v1,v2,v3,v4

!!!!! getting the values of interpolator variable
t1=VtoT(v1,this%v1Min,this%p11)
t2=VtoT2(v2,this%v2Min,this%p21)
t3=VtoT(v3,this%v3Min,this%p31)
t4=VtoT(v4,this%v4Min,this%p41)

!write(*,*) "got:",t1,t2,t3,t4
!write(*,*) "--->>",v2,this%v2Min,this%p21

!!!!! seting sector of the grid
if(t1<1) then
n1=0
elseif(t1>this%N1-1) then
n1=this%N1-2
else
n1=int(t1)-1
end if

if(t2<1) then
n2=0
elseif(t2>this%N2-1) then
n2=this%N2-2
else
n2=int(t2)-1
end if

if(t3<1) then
n3=0
elseif(t3>this%N3-1) then
n3=this%N3-2
else
n3=int(t3)-1
end if

if(t4<1) then
n4=0
elseif(t4>this%N4-1) then
n4=this%N4-2
else
n4=int(t4)-1
end if
!!!!!!!! basically the sector is n1,n1+1,n1+2,n1+3

!!!! extract the  subgrid
!!!! not that if the node is not used, it is computed
do i1=0,3
do i2=0,3
do i3=0,3
do i4=0,3
subGrid4D(i1,i2,i3,i4)=this%getValueAtNode(n1+i1,n2+i2,n3+i3,n4+i4)
end do
end do
end do
end do

!!!! interpolate
call interpolate4Dto3D(t1,n1,subGrid4D,subGrid3D)
call interpolate3Dto2D(t2,n2,subGrid3D,subGrid2D)
call interpolate2Dto1D(t3,n3,subGrid2D,subGrid1D)
call interpolate1DtoNUM(t4,n4,subGrid1D,getValue_this)


getValue_this=getValue_this/weigthFunction(v1,v2,v3,v4,this%process0)

end function getValue_this

!!!!!!! for the following subroutines:
!!!!! t = number on the grid
!!!!! n = number of the smallest node (i.e. n=0 for the lesser grid)
!!!!! p = is (0:3,....) subgid to interpolate

!!!!!! Neuvile algorithm for (:,:,:,:) grid 4D->3D
subroutine interpolate4Dto3D(t,n,p,res)
real(dp), dimension(0:3,0:3,0:3,0:3), intent(in)::p
real(dp), intent(out):: res(0:3,0:3,0:3)
real(dp), intent(in):: t
integer, intent(in)::n

real(dp):: p01(0:3,0:3,0:3),p12(0:3,0:3,0:3),p23(0:3,0:3,0:3), p012(0:3,0:3,0:3),p123(0:3,0:3,0:3)

p01(:,:,:)=(n+1-t)*p(0,:,:,:)+(t-n-0)*p(1,:,:,:)
p12(:,:,:)=(n+2-t)*p(1,:,:,:)+(t-n-1)*p(2,:,:,:)
p23(:,:,:)=(n+3-t)*p(2,:,:,:)+(t-n-2)*p(3,:,:,:)

p012=(n+2-t)*p01+(t-n-0)*p12
p123=(n+3-t)*p12+(t-n-1)*p23

res=((n+3-t)*p012+(t-n-0)*p123)/6 !!!!! division by 6 is because of 1/2 in the second step and 1/3 in the last step


end subroutine interpolate4Dto3D

!!!!!! Neuvile algorithm for (:,:,:) grid 3D->2D
subroutine interpolate3Dto2D(t,n,p,res)
real(dp), dimension(0:3,0:3,0:3), intent(in)::p
real(dp), intent(out):: res(0:3,0:3)
real(dp), intent(in):: t
integer, intent(in)::n

real(dp):: p01(0:3,0:3),p12(0:3,0:3),p23(0:3,0:3), p012(0:3,0:3),p123(0:3,0:3)

p01(:,:)=(n+1-t)*p(0,:,:)+(t-n-0)*p(1,:,:)
p12(:,:)=(n+2-t)*p(1,:,:)+(t-n-1)*p(2,:,:)
p23(:,:)=(n+3-t)*p(2,:,:)+(t-n-2)*p(3,:,:)

p012=(n+2-t)*p01+(t-n-0)*p12
p123=(n+3-t)*p12+(t-n-1)*p23

res=((n+3-t)*p012+(t-n-0)*p123)/6 !!!!! division by 6 is because of 1/2 in the second step and 1/3 in the last step

end subroutine interpolate3Dto2D

!!!!!! Neuvile algorithm for (:,:) grid 2D->1D
subroutine interpolate2Dto1D(t,n,p,res)
real(dp), dimension(0:3,0:3), intent(in)::p
real(dp), intent(out):: res(0:3)
real(dp), intent(in):: t
integer, intent(in)::n

real(dp):: p01(0:3),p12(0:3),p23(0:3), p012(0:3),p123(0:3)

p01(:)=(n+1-t)*p(0,:)+(t-n-0)*p(1,:)
p12(:)=(n+2-t)*p(1,:)+(t-n-1)*p(2,:)
p23(:)=(n+3-t)*p(2,:)+(t-n-2)*p(3,:)

p012=(n+2-t)*p01+(t-n-0)*p12
p123=(n+3-t)*p12+(t-n-1)*p23

res=((n+3-t)*p012+(t-n-0)*p123)/6 !!!!! division by 6 is because of 1/2 in the second step and 1/3 in the last step

end subroutine interpolate2Dto1D

!!!!!! Neuvile algorithm for (:) grid 1D->RESULT
subroutine interpolate1DtoNUM(t,n,p,res)
real(dp), dimension(0:3), intent(in)::p
real(dp), intent(out):: res
real(dp), intent(in):: t
integer, intent(in)::n

real(dp):: p01,p12,p23, p012,p123

p01=(n+1-t)*p(0)+(t-n-0)*p(1)
p12=(n+2-t)*p(1)+(t-n-1)*p(2)
p23=(n+3-t)*p(2)+(t-n-2)*p(3)

p012=(n+2-t)*p01+(t-n-0)*p12
p123=(n+3-t)*p12+(t-n-1)*p23

res=((n+3-t)*p012+(t-n-0)*p123)/6 !!!!! division by 6 is because of 1/2 in the second step and 1/3 in the last step

end subroutine interpolate1DtoNUM

end module aTMDe_xGrid
