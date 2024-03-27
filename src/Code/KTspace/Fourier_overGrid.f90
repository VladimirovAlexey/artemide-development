!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.0
!
!    This file contains the module, which is common for all TMD-evaluation modules
!
!                A.Vladimirov (16.03.2024)
!---------------------------------------------------------------------------------------

!!!!
!!!! This module performes the Fourier tranfomation to kT-space.
!!!! It does it by the Levin method, and transform the Chebyshev grid to the Chebyshev grid.
!!!! The result is the grid in the (x,kT,Q)-space [4D!]
!!!!

!     In the file that uses it add
!module NAME
! INCLUDE this_file
!end module NAME
!

!module Fourier_overGrid

use aTMDe_Numerics
use IO_functions
use InverseMatrix, only: Inverse

implicit none

private

character(len=8),parameter :: moduleName="tw2-grid"
character(len=11):: parentModuleName
integer::outputLevel


!!!!---------- Variables about the grid in b-space
!!!! bValues is the list of b-values required for the transform
real(dp),allocatable::bValues(:)

!!!!---------- Variables about the grid in k-space
!!!! kRanges are the list of values of subgrids for kT
!!!! Nodes are the list of values of nodes in the terms of T=[-1,1]
!!!! NodeFactors are the list of factors (-1)^i/2^%., there %=1 for i=0,NUM, and 0 otherwice
real(dp),allocatable::kRanges(:),kNodes(:),kNodeFactors(:)
!!!! number of Subgrids
integer::numKsubgrids
!!!! number of nodes, it is the same for all subgrids of given class
integer::kGridSize

!!!!---------- Variables about the grid in k-space
!!!! kRanges are the list of values of subgrids for kT
!!!! Nodes are the list of values of nodes in the terms of T=[-1,1]
!!!! NodeFactors are the list of factors (-1)^i/2^%., there %=1 for i=0,NUM, and 0 otherwice
real(dp),allocatable::xRanges(:),xNodes(:),xNodeFactors(:)
!!!! number of Subgrids
integer::numXsubgrids
!!!! number of nodes, it is the same for all subgrids of given class
integer::xGridSize


!!!!---------- General parameters
!!!! number of hadrons
integer::numH
!!!! include Gluon
logical::withGluon
!!!! utmost values of the grids
real(dp)::QMIN,QMAX,kMIN,kMAX,xMIN

real(dp)::zero=10.d-12
!!!! first zero of bessel function J0,J1,...,J5
real(dp),dimension(0:5),parameter::besselZERO=(/&
2.4048255576957727686216318793265_dp, 3.8317059702075123156144358863082_dp, &
5.1356223018406825563014016901378_dp, 6.3801618959239835062366146419427_dp, &
7.5883424345038043850696300079857_dp, 8.7714838159599540191228671334096_dp/)

!!!! Intervals are the lists of (u_{k+1}-u_k)/2 for subgrids
!!!! Means are the lists of (u_{k+1}+u_k)/2 for subgrids
!!!! these are used to speed-up transformation from and to grids (the list are in the trnasformed variables)
real(dp),allocatable::xIntervals(:),xMeans(:),kIntervals(:),kMeans(:),QIntervals(:),QMeans(:)

public:: PrepareTransformationMatrix

contains


!!!!!! the parameters for q-grid and k-grid are to be read from the file (section for KT-grid)
!!!!!! the parameters for x-grid and b-grid are to be read from the file (section for b-Grid)
subroutine Initialize_Fourier_overGrid(name,outLevel,numH_in,withGluon_in)
integer,intent(in)::outLevel
character(*),intent(in)::name
integer,intent(in)::numH_in
logical,intent(in)::withGluon_in

!!!!---------- Variables about the grid in b-space
real(dp),allocatable::bRanges(:)
integer::numBsubgrids,bGridSize

integer::i

parentModuleName=name
outputLevel=outLevel

!!!processing input parameters
numXsubgrids=2
numBsubgrids=4
numKsubgrids=3
allocate(xRanges(0:numXsubgrids))
allocate(bRanges(0:numBsubgrids))
allocate(kRanges(0:numKsubgrids))

xRanges=(/0.001d0,0.2d0,1.d0/)
bRanges=(/0.00001d0,0.01d0,0.2d0,2.d0, 25.d0/)
kRanges=(/0.001d0,2.d0, 20.d0,100.d0/)
xGridSize=16
bGridSize=16
kGridSize=16
xMIN=xRanges(0)
kMIN=kRanges(0)
kMAX=kRanges(numKsubgrids)

numH=numH_in

withGluon=withGluon_in

!!!!allocation of lists

allocate(xNodes(0:xGridSize))
allocate(kNodes(0:kGridSize))
allocate(xNodeFactors(0:xGridSize))
allocate(kNodeFactors(0:kGridSize))

allocate(xIntervals(1:numXsubgrids),xMeans(1:numXsubgrids))
allocate(kIntervals(1:numKsubgrids),kMeans(1:numKsubgrids))

!!!!filing the working variables

xIntervals=(log(xRanges(1:numXsubgrids))-log(xRanges(0:numXsubgrids-1)))/2._dp
kIntervals=(log(kRanges(1:numKsubgrids))-log(kRanges(0:numKsubgrids-1)))/2._dp
xMeans=(log(xRanges(1:numXsubgrids))+log(xRanges(0:numXsubgrids-1)))/2._dp
kMeans=(log(kRanges(1:numKsubgrids))+log(kRanges(0:numKsubgrids-1)))/2._dp

xNodeFactors=1._dp

do i=0,xGridSize
  xNodes(i)=cos(i*pi/xGridSize)

  if(i==0 .or. i==xGridSize) then
    xNodeFactors(i)=xNodeFactors(i)/2
  end if
  if(mod(i,2)==1) xNodeFactors(i)=-xNodeFactors(i)

end do

kNodeFactors=1._dp

do i=0,kGridSize
  kNodes(i)=cos(i*pi/kGridSize)

  if(i==0 .or. i==kGridSize) then
    kNodeFactors(i)=kNodeFactors(i)/2
  end if
  if(mod(i,2)==1) kNodeFactors(i)=-kNodeFactors(i)
end do

call PrepareTransformationMatrix(bGridSize,bRanges)

end subroutine Initialize_Fourier_overGrid


!!!!! This subroutine creates the transformation matrix from b to kT space
!!!!! it incapsulates all variables and functions about the b-space, and makes the table of bValues
subroutine PrepareTransformationMatrix(bGridSize,bRanges)
integer,intent(in)::bGridSize
real(dp),intent(in)::bRanges(0:)!(0:numBsubgrids)
integer::numBsubgrids
real(dp)::bIntervals(1:size(bRanges)-1),bMeans(1:size(bRanges)-1),tNodes(0:bGridSize),bNodeFactors(0:bGridSize)

real(dp),dimension(0:bGridSize)::CCweight!!!!! weights of CC-quadrature

integer::i,j,k
real(dp)::FourierAtQ(1:((size(bRanges)-1)*(bGridSize+1)))
real(dp)::TEST(1:((size(bRanges)-1)*(bGridSize+1)))
real(dp)::dummy
real(dp)::qTEST

numBsubgrids=size(bRanges)-1

!!!!! filling working variables for b-grid
bIntervals=(log(bRanges(1:numBsubgrids))-log(bRanges(0:numBsubgrids-1)))/2._dp
bMeans=(log(bRanges(1:numBsubgrids))+log(bRanges(0:numBsubgrids-1)))/2._dp

bNodeFactors=1._dp

do i=0,bGridSize
  tNodes(i)=cos(i*pi/bGridSize)

  if(i==0 .or. i==bGridSize) then
    bNodeFactors(i)=bNodeFactors(i)/2
  end if
  if(mod(i,2)==1) bNodeFactors(i)=-bNodeFactors(i)
end do

!!!!!! make the list of required values of b
allocate(bValues(1:(numBsubgrids*(bGridSize+1))))
k=1
do i=1,numBsubgrids
do j=0,bGridSize
    bValues(k)=BfromNode(i,j)
    k=k+1
end do
end do

!!!! precompute CCweights
call CCweight_compute()

k=1
do i=1,numBsubgrids
do j=0,bGridSize
    dummy=BfromNode(i,j)
    TEST(k)=dummy*Exp(-2.1_dp*dummy)*(1+dummy**2*6.3_dp+dummy**6*1.4_dp)
    k=k+1
end do
end do

do i=1,100
qTEST=i/100._dp
FourierAtQ=FourierArrayAtQ(qTEST)
write(*,'("{",F8.4,",",F16.12,"},")',advance="no") qTEST,dot_product(FourierAtQ,TEST)
end do
write(*,*) " "
write(*,*) "---- "

do i=1,30
qTEST=i
FourierAtQ=FourierArrayAtQ(qTEST)
write(*,'("{",F8.4,",",F24.18,"},")',advance="no") qTEST,qTEST**2*dot_product(FourierAtQ,TEST)
end do
write(*,*) " "
write(*,*) "---- "

do i=10,100
qTEST=i*3
FourierAtQ=FourierArrayAtQ(qTEST)
write(*,'("{",F8.4,",",F24.18,"},")',advance="no") qTEST,qTEST**2*dot_product(FourierAtQ,TEST)
end do
write(*,*) " "
write(*,*) "---- "
stop

!
! do i=0,2*bGridSize+1
! do j=0,2*bGridSize+1
!     dummy=ElementOfLevinMatrix(i,j,3,100.d0)
!     write(*,'(F9.4)',advance="no") dummy
! end do
! write(*,*) " "
! end do
! ! stop
!
! inv=InverseLevinOperator(1,21.2d0)
! write(*,*) inv
! inv=InverseLevinOperator(2,21.2d0)
! write(*,*) inv
! inv=InverseLevinOperator(3,21.2d0)
! write(*,*) inv
! inv=InverseLevinOperator(4,21.2d0)
! write(*,*) inv

stop

contains

!!!!! values of b computed from the nodes
!!!!! n=subgrid, k=node
function BfromNode(n,k)
integer,intent(in)::n,k
real(dp)::BfromNode
BfromNode=exp(bIntervals(n)*tNodes(k)+bMeans(n))
end function BfromNode

!!!!! returns the array of weights, which should be multiplied by array of function at bNodes to get the Integral
!!!!! \int_bMIN^bMAX J_0(b q) f(b) db
!!!!! for the grids with bMAX*q<besselZERO*a, CC-quadrature is used, otherwise Levin
!!!!! the parameter a is any reasonable number
!!!!! --> if q is too small, the Levin matrix became unstable
!!!!! --> if q is too large, the CC quadrature does not work
!!!!! empirically I found that 1/2 workds good.
function FourierArrayAtQ(q)
real(dp),intent(in)::q
real(dp),dimension(1:((size(bRanges)-1)*(bGridSize+1)))::FourierArrayAtQ

real(dp)::inv(0:bGridSize)
integer::i

!!!! if the integral is not oscilating, one better use CC-quadrature.
do i=1,numBsubgrids
    if(q*bRanges(i)<besselZERO(0)*2) then
        inv=InverseCCOperator(i,q)
    else
        inv=InverseLevinOperator(i,q)
    end if
    FourierArrayAtQ(1+(i-1)*(bGridSize+1):i*(bGridSize+1))=inv(0:bGridSize)
end do

!!!! This is the correction for the extremely small-b
FourierArrayAtQ(bGridSize+1)=FourierArrayAtQ(bGridSize+1)+bessel_j1(bRanges(0)*q)/q

end function FourierArrayAtQ

!!!! returns the derivative matrix for Chebyshev grid (see 2.19 in [2112.09703])
!!!! it does not contains the prefactor 1/b/([B-A]/2) recuired for proper derivative
pure function DerivativeMatrix(j,k)
integer,intent(in):: j,k
real(dp)::DerivativeMatrix
real(dp)::dummy
if(j==k) then
    if(j==0) then
        DerivativeMatrix=real(2*bGridSize**2+1,dp)/6
    else if (j==bGridSize) then
        DerivativeMatrix=-real(2*bGridSize**2+1,dp)/6
    else
        dummy=Cos(j*pi/bGridSize)
        DerivativeMatrix=-dummy/(2*(1-dummy**2))
    end if
else
    if(j==0 .or. j==bGridSize) then
        dummy=2._dp
    else
        dummy=1._dp
    end if
    if(k==0 .or. k==bGridSize) then
        dummy=dummy/2._dp
    end if
    if(mod(k+j,2)==1) then
        dummy=-dummy
    end if
    DerivativeMatrix=dummy/(Cos(j*pi/bGridSize)-Cos(k*pi/bGridSize))
end if
end function DerivativeMatrix

!!!! computes the element of the Levin matrix for the b-subgrid n, and q
!!!! matrix is (0:2*bGridSize+1)
!!!! the matrix has block form
!!!! |  -qI  D-I/b |
!!!! |   D    qI   |
!!!! where D is derivative matrix (divided by interval*b[i] because of derivative) and I is identiy matrix
!!!! each block is (bGridSize+1) x (bGridSize+1)
!!!! This matrix is for {h2,h1} and additionally trnasposed (because the matrix with {h1,h2} cannot be Gaussian Eliminated simply)
function ElementOfLevinMatrix(i,j,n,q)
real(dp)::ElementOfLevinMatrix
real(dp),intent(in)::q
integer,intent(in)::n
integer::i,j


! !!! |  -qI  D-I/b |
! !!! |   D    qI   |
! if(i<bGridSize+1 .and. j<bGridSize+1) then                      !!!!! upper-left block
!     if(i==j) then
!         ElementOfLevinMatrix=-q
!     else
!         ElementOfLevinMatrix=0._dp
!     end if
! else if (i<bGridSize+1 .and. j>bGridSize) then                !!!!! upper-right block
!     ElementOfLevinMatrix=DerivativeMatrix(i,j-bGridSize-1)/bIntervals(n)/BfromNode(n,i)
!     if(i==j-bGridSize-1) then
!         ElementOfLevinMatrix=ElementOfLevinMatrix-1/BfromNode(n,i)
!     end if
!
! else if (i>bGridSize .and. j<bGridSize+1) then                !!!!! lower-left block
!     ElementOfLevinMatrix=DerivativeMatrix(i-bGridSize-1,j)/bIntervals(n)/BfromNode(n,i-bGridSize-1)
! else
!     if(i==j) then
!         ElementOfLevinMatrix=q
!     else
!         ElementOfLevinMatrix=0._dp
!     end if
! end if

!!!! |   D         q I   |
!!!! |   -q I    D-I/b   |
if(i<bGridSize+1 .and. j<bGridSize+1) then                      !!!!! upper-left block
    ElementOfLevinMatrix=DerivativeMatrix(i,j)/bIntervals(n)/BfromNode(n,i)
else if (i<bGridSize+1 .and. j>bGridSize) then                !!!!! upper-right block
    if(i==j-bGridSize-1) then
        ElementOfLevinMatrix=q
    else
        ElementOfLevinMatrix=0._dp
    end if

else if (i>bGridSize .and. j<bGridSize+1) then                !!!!! lower-left block
    if(i-bGridSize-1==j) then
        ElementOfLevinMatrix=-q
    else
        ElementOfLevinMatrix=0._dp
    end if
else
    ElementOfLevinMatrix=DerivativeMatrix(i-bGridSize-1,j-bGridSize-1)/bIntervals(n)/BfromNode(n,i-bGridSize-1)
    if(i==j) then
        ElementOfLevinMatrix=ElementOfLevinMatrix-1/BfromNode(n,i-bGridSize-1)
    end if
end if

end function ElementOfLevinMatrix

!!!! finds the inverse of the Levin Matrix by Gaussian elimination
!!!! then builds the array A_{i}J_1(B q)-A_{i}J_1(A q)+A_{i}J_0(B q)-A_{i}J_0(A q)
!!!! which corresponds to the integral \int_A^B J_0(qb)
function InverseLevinOperator(n,q)
real(dp),intent(in)::q
integer,intent(in)::n
integer::NUM
real(dp),dimension(0:bGridSize)::InverseLevinOperator

real(dp),dimension(0:2*bGridSize+1,0:2*bGridSize+1)::MM!!!! matrix to inverse
real(dp),dimension(0:2*bGridSize+1,0:2*bGridSize+1)::invMM!!!! result of inverstion
integer::i,j,k
real(dp)::sum_em

NUM=2*bGridSize+1


!!!! make a matrix to invert
do i=0,NUM
do j=0,NUM
    MM(i,j)=ElementOfLevinMatrix(i,j,n,q)
end do
end do

!!!!! inversion is done by Crout (realisation by https://github.com/Beliavsky/Matrix_Inversion)
invMM=inverse(MM)

!!!!! constract the transformation array
!!!!! NOTE: that t->1 corresponds to the lower limit, and t->-1 corresponds to the upper limit

InverseLevinOperator=invMM(0,0:bGridSize)*bessel_j0(BfromNode(n,0)*q) &
                    -invMM(bGridSize,0:bGridSize)*bessel_j0(BfromNode(n,bGridSize)*q) &
                    +invMM(bGridSize+1,0:bGridSize)*bessel_j1(BfromNode(n,0)*q) &
                    -invMM(NUM,0:bGridSize)*bessel_j1(BfromNode(n,bGridSize)*q)
end function InverseLevinOperator


!!!!! in some cases (very small q), the Levin-integration is not good.
!!!!! in these caes it is convenient to use Clenshaw Curtis since it uses the same nodes
!!!!! this function computes the cuadrature for the integral int J_0(bq)f[b] in the range n.
function InverseCCOperator(n,q)
real(dp),intent(in)::q
integer,intent(in)::n
real(dp),dimension(0:bGridSize)::InverseCCOperator

integer::i

do i=0,bGridSize
    InverseCCOperator(i)=bIntervals(n)*BfromNode(n,i)*bessel_j0(q*BfromNode(n,i))*CCweight(i)
end do

end function InverseCCOperator

!!!! precomputes the CCweights
!!!! the formula is
!!!! w_j=beta(j)4/N sum_{k=0,2,..}^n beta[k] Cos[j k pi/n]/(1-k^2)
!!!! where n is number of nodes, beta=1/2 for 0 or n.
subroutine CCweight_compute()
integer::j,k
do j=0,bGridSize
    CCweight(j)=0._dp
    do k=0,bGridSize,2
        if(k==0 .or. k==bGridSize) then
            CCweight(j)=CCweight(j)+Cos(j*k*pi/bGridSize)/2/(1._dp-k**2)
        else
            CCweight(j)=CCweight(j)+Cos(j*k*pi/bGridSize)/(1._dp-k**2)
        end if
    end do
    if(j==0 .or. j==bGridSize) CCweight(j)=CCweight(j)/2
    CCweight(j)=CCweight(j)*4._dp/bGridSize
end do
end subroutine CCweight_compute


end subroutine PrepareTransformationMatrix



!end module Fourier_overGrid
