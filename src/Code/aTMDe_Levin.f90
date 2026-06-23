!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.04
!
!    Contains the module that realizes the Fourier integration with Levin algorithms over a grid
!
!    This module performs the Fourier transformation to kT-space.
!    It does it by the Levin method, and transforms the Chebyshev grid to the Chebyshev grid.
!
!    It makes the Fourier of the function F with the interface (-5:5)
!    The transformations are defined as
!    for TMDtypeN=0:  \int_0^\infty db/(2pi) b J_0(b q) F(b)
!    for TMDtypeN=1:  M^2/kT\int_0^\infty db/(2pi) b^2 J_1(b q) F(b)
!
!    Has two functions:
!    1) returns value of Fourier at a specific point
!    2) Returns a grid in kT-space, which can be used by aTMDe_ktGrid)
!
!    v1. implementation !   (16.03.2024)
!    v2. Update to object structure (23.06.2026)
!                A.Vladimirov (23.06.2026)
!---------------------------------------------------------------------------------------

module aTMDe_Levin
use aTMDe_Numerics
use aTMDe_interfaces
use aTMDe_IO
use aTMDe_invMatrix
implicit none

private

character(len=5),parameter :: moduleName="Levin"

!!!! the first zero of bessel function J0,J1,...,J5
real(dp),dimension(0:5),parameter::besselZERO=(/&
2.4048255576957727686216318793265_dp, 3.8317059702075123156144358863082_dp, &
5.1356223018406825563014016901378_dp, 6.3801618959239835062366146419427_dp, &
7.5883424345038043850696300079857_dp, 8.7714838159599540191228671334096_dp/)

!!!!!! The object of Levin transform --------
type, public :: LevinIntegrator
    private

    character(len=11):: parentName
    integer::outputLevel

    !!!! General parameters
    !! mass parameter used as mass-scale
    real(dp)::TMDmass=1._dp
    !!!!! this is the order of Bessel-transform (IT IS STRICT FOR TMD; Must be 0, 1 in the present implementation)
    integer::TMDtypeN
    !!!!!
    logical:: prepareGrid=.false.

    !!!!---------- Variables about the grid in b-space
    !!!! bValues is the list of b-values from all subgrids, collected in the single list
    !!!! It is used for the transform, and has the shape (1:lengthOfbValues)
    real(dp),allocatable::bValues(:)
    !!!! the length of bValues is =(numBsubgrids*(bGridSize+1))
    integer::lengthOfbValues
    !!!! bRanges are the list of values of subgrids for b
    !!!! Nodes are the list of values of nodes in the terms of T=[-1,1]
    !!!! NodeFactors are the list of factors (-1)^i/2^%., where %=1 for i=0,NUM, and 0 otherwise
    real(dp),allocatable::bRanges(:),bNodes(:),bNodeFactors(:)
    !!!! number of Subgrids
    integer::numBsubgrids
    !!!! number of nodes, it is the same for all subgrids of given class
    integer::bGridSize
    !!!! Intervals are the lists of (u_{k+1}-u_k)/2 for subgrids
    !!!! Means are the lists of (u_{k+1}+u_k)/2 for subgrids
    !!!! these are used to speed up transformation from and to grids (the lists are in the transformed variables)
    real(dp),allocatable::bIntervals(:),bMeans(:)


    !!!!---------- Variables about the grid in k-space
    !!!! kRanges are the list of values of subgrids for kT
    !!!! Nodes are the list of values of nodes in the terms of T=[-1,1]
    !!!! NodeFactors are the list of factors (-1)^i/2^%., where %=1 for i=0,NUM, and 0 otherwise
    real(dp),allocatable::kRanges(:),kNodes(:),kNodeFactors(:)
    !!!! number of Subgrids
    integer::numKsubgrids
    !!!! number of nodes, it is the same for all subgrids of given class
    integer::kGridSize
    !!!! Intervals are the lists of (u_{k+1}-u_k)/2 for subgrids
    !!!! Means are the lists of (u_{k+1}+u_k)/2 for subgrids
    !!!! these are used to speed up transformation from and to grids (the lists are in the transformed variables)
    real(dp),allocatable::kIntervals(:),kMeans(:)

    !!!!--------- Variables concerning transformation to the kT-space
    !!!! The matrix or arrays (for each node in kT) to transform f(b)
    real(dp),allocatable::TransformationArray(:,:,:)

    !!!!! weights of CC-quadrature
    real(dp),allocatable::CCweight(:)

contains
  procedure,public::Fourier_atPoint =>Fourier_Levin
  procedure,public::Fourier_array  =>Fourier_Levin_array
end type

interface LevinIntegrator
    procedure :: constructor
end interface LevinIntegrator

contains


!!!!!! Initialization of parameters from the const-file
!!! path: path to INI-file
!!! moduleLine: "*12  " the line for the module input
!!! gridLine: "*F   " the line for the KT-grid input within the module section
!!! name: name of the parent module (for messages)
!!! outLevel: output level (for messages)
!!! TMDt: is the type of KT-transform for TMD (n=0, b J_0;  n=1 b^2 J_1)
function constructor(path,moduleLine,gridLine,name,outLevel,TMDt) result(this)
type(LevinIntegrator)::this
character(len=*),intent(in)::path
character(len=5),intent(in)::moduleLine,gridLine
integer,intent(in)::outLevel,TMDt
character(*),intent(in)::name

integer::i,j,k

this%parentName=name
this%outputLevel=outLevel

this%TMDtypeN=TMDt

if(this%TMDtypeN/=0 .and. this%TMDtypeN/=1) then
    ERROR STOP ErrorString(&
    "TMDtype for Levin transform is defined only for 0 or 1, while called "//int4ToStr(this%TMDtypeN),this%parentName,moduleName)
end if

!!!! read input about b and kT-spaces
OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    !------------- General parameters
    call MoveTO(51,'*0   ')
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p2  ')
    read(51,*) this%TMDmass
    !-------------Parameters of grid in kT
    call MoveTO(51,moduleLine)
    call MoveTO(51,gridLine)
    call MoveTO(51,'*p1  ')
    read(51,*) this%prepareGrid
    call MoveTO(51,'*p6  ')
    read(51,*) this%numKsubgrids
    allocate(this%kRanges(0:this%numKsubgrids))
    call MoveTO(51,'*p7  ')
    read(51,*) this%kRanges
    call MoveTO(51,'*p8  ')
    read(51,*) this%kGridSize

    !-------------Parameters of grid in b
    call MoveTO(51,'*p12 ')
    read(51,*) this%numBsubgrids
    allocate(this%bRanges(0:this%numBsubgrids))
    call MoveTO(51,'*p13 ')
    read(51,*) this%bRanges
    call MoveTO(51,'*p14 ')
    read(51,*) this%bGridSize
CLOSE (51, STATUS='KEEP')


!!!!allocation of lists
allocate(this%kNodes(0:this%kGridSize))
allocate(this%kNodeFactors(0:this%kGridSize))
allocate(this%kIntervals(1:this%numKsubgrids),this%kMeans(1:this%numKsubgrids))

allocate(this%bNodes(0:this%bGridSize))
allocate(this%bNodeFactors(0:this%bGridSize))
allocate(this%bIntervals(1:this%numBsubgrids),this%bMeans(1:this%numBsubgrids))

allocate(this%CCweight(0:this%bGridSize))

!    It precomputes the variable bValues (which is a list of b's)
!    And the TransformationArray. Applying transformation array to the f(bValues), one gets a grid in kT-space

!!!!filing the working variables for K-grid
this%kIntervals=(log(this%kRanges(1:this%numKsubgrids))-log(this%kRanges(0:this%numKsubgrids-1)))/2._dp
this%kMeans=(log(this%kRanges(1:this%numKsubgrids))+log(this%kRanges(0:this%numKsubgrids-1)))/2._dp

this%kNodeFactors=1._dp

do i=0,this%kGridSize
  this%kNodes(i)=cos(i*pi/this%kGridSize)

  if(i==0 .or. i==this%kGridSize) then
    this%kNodeFactors(i)=this%kNodeFactors(i)/2
  end if
  if(mod(i,2)==1) this%kNodeFactors(i)=-this%kNodeFactors(i)
end do

!!!!! filling working variables for b-grid
this%bIntervals=(log(this%bRanges(1:this%numBsubgrids))-log(this%bRanges(0:this%numBsubgrids-1)))/2._dp
this%bMeans=(log(this%bRanges(1:this%numBsubgrids))+log(this%bRanges(0:this%numBsubgrids-1)))/2._dp

this%bNodeFactors=1._dp

do i=0,this%bGridSize
  this%bNodes(i)=cos(i*pi/this%bGridSize)

  if(i==0 .or. i==this%bGridSize) then
    this%bNodeFactors(i)=this%bNodeFactors(i)/2
  end if
  if(mod(i,2)==1) this%bNodeFactors(i)=-this%bNodeFactors(i)
end do

!!!!!! make the list of required values of b
this%lengthOfbValues=this%numBsubgrids*(this%bGridSize+1)
allocate(this%bValues(1:this%lengthOfbValues))
!!! precompute bValues
k=1
do i=1,this%numBsubgrids
do j=0,this%bGridSize
    this%bValues(k)=BfromNode(this,i,j)
    k=k+1
end do
end do

!!!! precompute CCweights
call CCweight_compute(this)


if(this%prepareGrid) then
    allocate(this%TransformationArray(1:this%numKsubgrids,0:this%kGridSize,1:this%lengthOfbValues))
    call PrepareTransformationMatrix(this)
    if(this%outputLevel>=2) &
     write(*,*) "Initialization of tables for Hankel transform over grid is done for "//trim(this%parentName)
end if

end function constructor

!!!!! values of K computed from the nodes
!!!!! n=subgrid, k=node
pure function KfromNode(this,n,k)
type(LevinIntegrator),intent(in)::this
integer,intent(in)::n,k
real(dp)::KfromNode
KfromNode=exp(this%KIntervals(n)*this%KNodes(k)+this%KMeans(n))
end function KfromNode

!!!!! values of b computed from the nodes
!!!!! n=subgrid, k=node
pure function BfromNode(this,n,k)
type(LevinIntegrator),intent(in)::this
integer,intent(in)::n,k
real(dp)::BfromNode
BfromNode=exp(this%bIntervals(n)*this%bNodes(k)+this%bMeans(n))
end function BfromNode

!!!!! This subroutine creates the transformation matrix from b to kT space
!!!!! the transformation is
!!!!! for TMDtypeN=0:  kT^2\int_0^\infty db/(2pi) J_0(b kT) [b F(b)]
!!!!! for TMDtypeN=1:  kT^2 M^2/kT\int_0^\infty db/(2pi) J_1(b kT) [b^2 F(b)]
!!!!! It works over the grid
subroutine PrepareTransformationMatrix(this)
type(LevinIntegrator),intent(inout)::this
integer::i,j

    !!!!! actually computation of the arrays for Fourier transform
SELECT CASE(this%TMDtypeN)
    CASE(0) !!!! uTMDPDF, and similar
        do i=1,this%numKsubgrids
        do j=0,this%kGridSize
            this%TransformationArray(i,j,:)=FourierArrayAtQ(this,KfromNode(this,i,j))/pix2*KfromNode(this,i,j)**2
        end do
        end do
    CASE(1) !!!! SiversTMDPDF, and similar
        do i=1,this%numKsubgrids
        do j=0,this%kGridSize
            this%TransformationArray(i,j,:)=FourierArrayAtQ(this,KfromNode(this,i,j))/pix2*(this%TMDMass**2)*KfromNode(this,i,j)
        end do
        end do
    CASE DEFAULT
        ERROR STOP ErrorString("Fourier_Levin: Unknown TMDtype. Presently implemeted only types 0,1",this%parentName,moduleName)
END SELECT

end subroutine PrepareTransformationMatrix

!!!!!! make the Fourier of the function F with the interface (-5:5)
!!!!!! at the point kT (in GeV)
!!!!!! the transformation is
!!!!!! for TMDtypeN=0:  \int_0^\infty db/(2pi) b J_0(b q) F(b)
!!!!!! for TMDtypeN=1:  M^2/kT\int_0^\infty db/(2pi) b^2 J_1(b q) F(b)
function Fourier_Levin(this,F,kT)
class(LevinIntegrator),intent(in)::this
real(dp),intent(in)::kT
real(dp),dimension(-5:5)::Fourier_Levin
procedure(func_1D_array5)::F
integer::i
real(dp),dimension(1:this%lengthOfbValues)::Levin_Array

Levin_Array=FourierArrayAtQ(this,kT)

SELECT CASE(this%TMDtypeN)
    CASE(0) !!! uTMDPDF, and similar
        Fourier_Levin=0._dp
        do i=1,this%lengthOfbValues
        Fourier_Levin=Fourier_Levin+Levin_Array(i)*F(this%bValues(i))*this%bValues(i)
        end do
        Fourier_Levin=Fourier_Levin/pix2

    CASE(1) !!! SiversTMDPDF, and similar
        Fourier_Levin=0._dp
        do i=1,this%lengthOfbValues
        Fourier_Levin=Fourier_Levin+Levin_Array(i)*F(this%bValues(i))*this%bValues(i)*this%bValues(i)
        end do
        Fourier_Levin=Fourier_Levin/pix2*this%TMDMass*this%TMDMass/kT
    CASE DEFAULT
        ERROR STOP ErrorString("Fourier_Levin: Unknown TMDtype. Presently impleneted only types 0,1",this%parentName,moduleName)

END SELECT

end function Fourier_Levin

!!!!!! make the Fourier of the function F with the interface (-5:5)
!!!!!! at the modes of kT-grid. Returns the list (1:numKsubgrids, 0:kGridSize,-5:5)
!!!!!! the transformation is
!!!!!! TMDtypeN=0:   kT^2\int_0^\infty db/(2pi) b J_0(b kT) F(b)
!!!!!! TMDtypeN=1:   kT^2 M**2/kT\int_0^\infty db/(2pi) b^2 J_1(b kT) F(b)
!!!!!! IMPORTANT THE RESULT IS MULTIPLIED by kT^2
!!!!!! Requires prepareGrid=.true. in the INI file; otherwise TransformationArray is not allocated.
function Fourier_Levin_array(this,F2)
class(LevinIntegrator),intent(in)::this
real(dp),dimension(1:this%numKsubgrids, 0:this%kGridSize,-5:5)::Fourier_Levin_array
procedure(func_1D_array5)::F2
integer::i,j,ff
real(dp),dimension(1:this%lengthOfbValues,-5:5)::Function_Array

!!!! request array
if(this%TMDtypeN==0) then
    do i=1,this%lengthOfbValues
    Function_Array(i,:)=F2(this%bValues(i))*this%bValues(i)
    end do
else if(this%TMDtypeN==1) then
    do i=1,this%lengthOfbValues
    Function_Array(i,:)=F2(this%bValues(i))*this%bValues(i)*this%bValues(i)
    end do
!!!!!
end if

!!!!! TMDtypeN=0: factor kT^2/2pi is taken into account in the TransformationArray
!!!!! TMDtypeN=1: factor kT^2/2pi M^2/kT is taken into account in the TransformationArray
do i=1,this%numKsubgrids
do j=0,this%kGridSize
do ff=-5,5
    Fourier_Levin_array(i,j,ff)=dot_product(this%TransformationArray(i,j,:),Function_Array(:,ff))
end do
end do
end do

end function Fourier_Levin_array

!!!!! returns the array of weights, which should be multiplied by array of function at bNodes to get the Integral
!!!!! for TMDtypeN=0: \int_bMIN^bMAX J_0(b q) f(b) db
!!!!! for TMDtypeN=1: \int_bMIN^bMAX J_1(b q) f(b) db
!!!!! for the grids with bMAX*q<besselZERO*a, CC-quadrature is used, otherwise Levin
!!!!! the parameter a is any reasonable number, which is set to a=bGridSize**2/128
!!!!! --> if q is too small, the Levin matrix becomes unstable
!!!!! --> if q is too large, the CC quadrature does not work
!!!!! empirically I found that 1/2 works well.
function FourierArrayAtQ(this,q)
type(LevinIntegrator),intent(in)::this
real(dp),intent(in)::q
real(dp),dimension(1:this%lengthOfbValues)::FourierArrayAtQ

real(dp)::inv(0:this%bGridSize)
integer::i
real(dp)::valueTOswitch

!!!! this is value at which I change CC to Levin.
!!!! if (bq) is very small, the integral is not oscillating, and is better estimated by CC
!!!! It should be compared with BesselZero.
!!!! Empirically, I found that the best threshold is:
valueTOswitch=besselZERO(this%TMDtypeN)*this%bGridSize**2/128.d0

!!!! if the integral is not oscillating, one should better use CC-quadrature.
do i=1,this%numBsubgrids
    if(q*this%bRanges(i)<valueTOswitch) then
        inv=InverseCCOperator(i,q)
    else
        inv=InverseLevinOperator(i,q)
    end if
    FourierArrayAtQ(1+(i-1)*(this%bGridSize+1):i*(this%bGridSize+1))=inv(0:this%bGridSize)
end do

!!!! This is the correction for the extremely small-b
!!!! In the case of type=0, this adds the missing [0, bRanges(0)] contribution, using the analytic result int_0^A b J_0(qb) db = A J_1(qA)/q.
!!!! For type=1 (and higher) integrals such correction is not needed, because the integrand b^2 J_1(qb) -> 0 as b->0
if(this%TMDtypeN==0) then
    FourierArrayAtQ(this%bGridSize+1)=FourierArrayAtQ(this%bGridSize+1)+bessel_j1(this%bRanges(0)*q)/q
end if

contains

!!!! returns the derivative matrix for Chebyshev grid (see 2.19 in [2112.09703])
!!!! it does not contain the prefactor 1/b/([B-A]/2) required for proper derivative
pure function DerivativeMatrix(j,k)
integer,intent(in):: j,k
real(dp)::DerivativeMatrix
real(dp)::dummy
if(j==k) then
    if(j==0) then
        DerivativeMatrix=real(2*this%bGridSize**2+1,dp)/6
    else if (j==this%bGridSize) then
        DerivativeMatrix=-real(2*this%bGridSize**2+1,dp)/6
    else
        dummy=Cos(j*pi/this%bGridSize)
        DerivativeMatrix=-dummy/(2*(1-dummy**2))
    end if
else
    if(j==0 .or. j==this%bGridSize) then
        dummy=2._dp
    else
        dummy=1._dp
    end if
    if(k==0 .or. k==this%bGridSize) then
        dummy=dummy/2._dp
    end if
    if(mod(k+j,2)==1) then
        dummy=-dummy
    end if
    DerivativeMatrix=dummy/(Cos(j*pi/this%bGridSize)-Cos(k*pi/this%bGridSize))
end if
end function DerivativeMatrix

!!!! computes the element of the Levin matrix for the b-subgrid n, and q
!!!! matrix is (0:2*bGridSize+1)
!!!! the matrix has block form
!!!! |   D         q I   |
!!!! |   -q I    D-I/b   |
!!!! where D is derivative matrix (divided by interval*b[i] because of derivative) and I is identity matrix
!!!! each block is (bGridSize+1) x (bGridSize+1)
!!!! This matrix is for {h2,h1} and additionally transposed (because the matrix with {h1,h2} cannot be Gaussian Eliminated simply)
function ElementOfLevinMatrix(i,j,n,q)
real(dp)::ElementOfLevinMatrix
real(dp),intent(in)::q
integer,intent(in)::n
integer::i,j

!!!! |   D         q I   |
!!!! |   -q I    D-I/b   |
if(i<this%bGridSize+1 .and. j<this%bGridSize+1) then                      !!!!! upper-left block
    ElementOfLevinMatrix=DerivativeMatrix(i,j)/this%bIntervals(n)/BfromNode(this,n,i)
else if (i<this%bGridSize+1 .and. j>this%bGridSize) then                !!!!! upper-right block
    if(i==j-this%bGridSize-1) then
        ElementOfLevinMatrix=q
    else
        ElementOfLevinMatrix=0._dp
    end if

else if (i>this%bGridSize .and. j<this%bGridSize+1) then                !!!!! lower-left block
    if(i-this%bGridSize-1==j) then
        ElementOfLevinMatrix=-q
    else
        ElementOfLevinMatrix=0._dp
    end if
else
    ElementOfLevinMatrix=&
        DerivativeMatrix(i-this%bGridSize-1,j-this%bGridSize-1)/this%bIntervals(n)/BfromNode(this,n,i-this%bGridSize-1)
    if(i==j) then
        ElementOfLevinMatrix=ElementOfLevinMatrix-1/BfromNode(this,n,i-this%bGridSize-1)
    end if
end if

end function ElementOfLevinMatrix

!!!! finds the inverse of the Levin Matrix
!!!! then builds the array A_{i}J_1(B q)-A_{i}J_1(A q)+A_{i}J_0(B q)-A_{i}J_0(A q)
!!!! which corresponds to the integral \int_A^B J_0(qb) f(b) db
function InverseLevinOperator(n,q)
real(dp),intent(in)::q
integer,intent(in)::n
integer::NUM
real(dp),dimension(0:this%bGridSize)::InverseLevinOperator

real(dp),dimension(0:2*this%bGridSize+1,0:2*this%bGridSize+1)::MM!!!! matrix to inverse
real(dp),dimension(0:2*this%bGridSize+1,0:2*this%bGridSize+1)::invMM!!!! its inverstion
integer::i,j

NUM=2*this%bGridSize+1


!!!! make a matrix to invert
do i=0,NUM
do j=0,NUM
    MM(i,j)=ElementOfLevinMatrix(i,j,n,q)
end do
end do

!!!!! inversion is done by Crout (realisation by https://github.com/Beliavsky/Matrix_Inversion)
invMM=inverse(MM)

!!!!! construct the transformation array
!!!!! NOTE: that t->1 corresponds to the upper limit (cos(0)=1), and t->-1 corresponds to the lower limit (cos(pi)=-1)

!!!! for type=0, I take the upper part of the matrix
!!!! for type=1, I take the lower part of the matrix
if(this%TMDtypeN==0) then
    InverseLevinOperator=invMM(0,0:this%bGridSize)*bessel_j0(BfromNode(this,n,0)*q) &
                        -invMM(this%bGridSize,0:this%bGridSize)*bessel_j0(BfromNode(this,n,this%bGridSize)*q) &
                        +invMM(this%bGridSize+1,0:this%bGridSize)*bessel_j1(BfromNode(this,n,0)*q) &
                        -invMM(NUM,0:this%bGridSize)*bessel_j1(BfromNode(this,n,this%bGridSize)*q)
else if(this%TMDtypeN==1) then
    InverseLevinOperator=invMM(0,this%bGridSize+1:NUM)*bessel_j0(BfromNode(this,n,0)*q) &
                        -invMM(this%bGridSize,this%bGridSize+1:NUM)*bessel_j0(BfromNode(this,n,this%bGridSize)*q) &
                        +invMM(this%bGridSize+1,this%bGridSize+1:NUM)*bessel_j1(BfromNode(this,n,0)*q) &
                        -invMM(NUM,this%bGridSize+1:NUM)*bessel_j1(BfromNode(this,n,this%bGridSize)*q)
!!!
end if
end function InverseLevinOperator


!!!!! in some cases (very small q), the Levin-integration is not good.
!!!!! in these cases it is convenient to use Clenshaw Curtis since it uses the same nodes
!!!!! this function computes the quadrature for the integral int J_0,1(bq)f[b] in the range n.
function InverseCCOperator(n,q)
real(dp),intent(in)::q
integer,intent(in)::n
real(dp),dimension(0:this%bGridSize)::InverseCCOperator

integer::i

if(this%TMDtypeN==0) then
    do i=0,this%bGridSize
        InverseCCOperator(i)=this%bIntervals(n)*BfromNode(this,n,i)*bessel_j0(q*BfromNode(this,n,i))*this%CCweight(i)
    end do
else if(this%TMDtypeN==1) then
    do i=0,this%bGridSize
        InverseCCOperator(i)=this%bIntervals(n)*BfromNode(this,n,i)*bessel_j1(q*BfromNode(this,n,i))*this%CCweight(i)
    end do
!!!!!
end if

end function InverseCCOperator

end function FourierArrayAtQ


!!!! precomputes the CCweights
!!!! the formula is
!!!! w_j=beta(j)*4/N sum_{k=0,2,..}^n beta[k] Cos[j k pi/n]/(1-k^2)
!!!! where n is number of nodes, beta=1/2 for 0 or n.
subroutine CCweight_compute(this)
type(LevinIntegrator),intent(inout)::this
integer::j,k
do j=0,this%bGridSize
    this%CCweight(j)=0._dp
    do k=0,this%bGridSize,2
        if(k==0 .or. k==this%bGridSize) then
            this%CCweight(j)=this%CCweight(j)+Cos(j*k*pi/this%bGridSize)/2/(1._dp-k**2)
        else
            this%CCweight(j)=this%CCweight(j)+Cos(j*k*pi/this%bGridSize)/(1._dp-k**2)
        end if
    end do
    if(j==0 .or. j==this%bGridSize) this%CCweight(j)=this%CCweight(j)/2
    this%CCweight(j)=this%CCweight(j)*4._dp/this%bGridSize
end do
end subroutine CCweight_compute

end module aTMDe_Levin
