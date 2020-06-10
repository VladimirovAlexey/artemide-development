!!!---------------------------------------------------------------------
!!!
!!! The module that creates and reads the spline grids
!!!
!!!
!!!---------------------------------------------------------------------

module SplineGrid

implicit none

public

contains

subroutine CreateGrid1D(f1D,S,xMin,xMax,delta)
    real*8::f1D,x
    real*8,intent(in)::xMin,xMax !!! grid boundaries
    real*8,intent(in)::delta    !!! relative precision (fExact(x)-fGrid(x))/f(x)<delta
    real*8,allocatable,intent(out)::S(:,:)
    
    real*8,allocatable::xL(:),yL(:)
    integer::i,nCurrent,j
    integer,parameter::MaxIteration=100
    logical::toUpdateAgain
    
    !!!! the initial grid is made by equidistant lattice
    nCurrent=10
    allocate(xL(0:nCurrent),yL(0:nCurrent))
    
    !!!create nodes
    do i=0,nCurrent
        xL(i)=xMin*(1d0-i/real(nCurrent))+xMax*i/real(nCurrent)
        yL(i)=f1D(xL(i))
    end do
    
    !!!initial grid
    call CreateGrid1D_byNodes(S,xL,yL)
    
    !!!! iterate grid
    do j=1,MaxIteration
        deallocate(xL,yL)
        !!! intial grid complete
        
        nCurrent=size(S,1)-1
        !!! create intermidiate points
        allocate(xL(0:nCurrent-1),yL(0:nCurrent-1))
        do i=0,nCurrent-1
            xL(i)=(S(i,1)+S(i+1,1))/2d0
            yL(i)=f1D(xL(i))
        end do
        
        call TestAndIterateGrid(S,xL,yL,delta,toUpdateAgain)
        
        write(*,*) nCurrent,"-->",size(S,1)-1
        
        if(.not.toUpdateAgain) exit
    end do
    deallocate(xL,yL)
        
        
    write(*,*) SplineInterpolate1D(S,3.5d0)
    write(*,*) SplineInterpolate1D(S,4.6d0)
    write(*,*) SplineInterpolate1D(S,6.7d0)

end subroutine CreateGrid1D

!!!!! create the spline using the given nodes 
!!!!! the algorithm is taken from the text book
!!!!! S is output spline (...,(xi,ai,bi,ci,di),...), 
!!!!! the last spline in S is the fake (xn,yn,0,0,0), it is needed only to remember xn
!!!!! xArray,yArray is the lists of (x,f(x))
subroutine CreateGrid1D_byNodes(S,xArray,yArray)
    real*8,intent(in)::xArray(0:),yArray(0:)
    real*8,allocatable,intent(out)::S(:,:)
    real*8,allocatable::a(:),b(:)
    integer::n,i
    
    n=size(xArray)-1
    
    allocate(a(0:n),b(0:n))
    
    a(0:n)=yArray(0:n)
    b(0:n-1)=(yArray(1:n)-yArray(0:n-1))/(xArray(1:n)-xArray(0:n-1))
!     do i=0,n-1
!         b(i)=(yArray(i+1)-yArray(i))/(xArray(i+1)-xArray(i))
!     end do
!     b(n)=0
    
    allocate(S(0:n,1:3))
    do i=0,n-1
        S(i,1:3)=(/xArray(i),a(i),b(i)/)
    end do
    S(n,1:3)=(/xArray(n),a(n),0d0/)
    deallocate(a,b)
end subroutine CreateGrid1D_byNodes

!!!!! Return the interpolated value from the spline set
!!!!! for the values outside the grid extrapolate by the last constant
function SplineInterpolate1D(S,x)
    real*8::SplineInterpolate1D
    real*8,intent(in)::S(0:,1:)
    real*8,intent(in)::x
    integer::n,m,i,j,k
    real*8::xx
    
    n=size(S,1)-1
    
    !!! for the values outside the boundary of grid use the constants
    if(x>S(n,1)) then
        SplineInterpolate1D=S(n,2)
        return
    end if
    if(x<S(0,1)) then
        SplineInterpolate1D=S(0,2)
        return
    end if
    
    !!! binary search for for m, such that x(m) <= x <= x(m+1)
    i = 0
    j = n
    do while (j > i+1)
        k = (i+j)/2
        if(x < S(k,1)) then
            j=k
        else
            i=k
        end if
    end do
    
    xx=x-S(i,1)
    
    SplineInterpolate1D=S(i,2)+xx*S(i,3)

end function SplineInterpolate1D

!!!!! The routine which updates the grid
!!!!! it takes the input grid, and suggested values of new points
!!!!! It checks whenever the function on new points far from the original
!!!!! if yes the new list of nodes is made and the values of grid S updated.
!!!!! The logical argument IsUpdated is set to .true. 
subroutine TestAndIterateGrid(S,xL,yL,delta,IsUpdated)
    real*8,allocatable,intent(inout)::S(:,:)
    real*8,intent(in)::xL(0:),yL(0:),delta
    logical::IsUpdated
    
    integer::n,i,j
    logical,allocatable::updatIt(:)
    real*8,allocatable::xNew(:),yNew(:)
    real*8::ff
    
    
    n=size(S,1)-1
    
!     write(*,*) '-------------ITER------   n=',n
    !!! check which segments requare the update
    allocate(updatIt(0:n-1))
    j=0
    do i=0,n-1
        ff=SplineInterpolate1D(S,xL(i))
!         write(*,*) xL(i),yL(i)-ff,delta,Abs(yL(i)-ff)<delta
        if(Abs(yL(i)-ff)<delta) then
            updatIt(i)=.false.
        else
            updatIt(i)=.true.
            j=j+1
        end if
    end do
!     write(*,*) 'TO UPDATE',j,updatIt
    
    !!! if there is no update, return
    if(j==0) then
        IsUpdated=.false.
        return
    end if
    
    IsUpdated=.true.
    !!! create new x and y's
    allocate(xNew(0:n+j),yNew(0:n+j))
    j=0
    do i=0,n-1
        xNew(j)=S(i,1)
        yNew(j)=S(i,2)
        j=j+1
        if(updatIt(i)) then
            xNew(j)=xL(i)
            yNew(j)=yL(i)
            j=j+1
        end if
    end do
    xNew(j)=S(n,1)
    yNew(j)=S(n,2)
    deallocate(updatIt)
    
!     write(*,*) '--',xNew
    
    !!! delete old grid, and create new
    deallocate(S)
    call CreateGrid1D_byNodes(S,xNew,yNew)
    
    deallocate(xNew,yNew)

end subroutine TestAndIterateGrid

end module SplineGrid

program test
use SplineGrid

real*8,allocatable::Sg(:,:)

call CreateGrid1D(f1,Sg,0d0,10d0,1d-4)

contains

 function f1(x)
 real*8::f1,x
  f1=0.33d0*x+(1d0+sin(2d0*x+3d0))/(x**2+0.54**2)+4d0*(1d0-0.05d0*x)*cos(sqrt(x+1d0))
 
 end function f1

end program test
