!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.00
!
!	Evaluation of the small-b OPE for uTMDPDF
!	
!	if you use this module please, quote ????.????
!
!	ver 3.00: release (AV, 18.07.2023)
!
!				A.Vladimirov (18.07.2023)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!---- Concept---
! This module
! * computes the tw2-convolution C x f[x,b] at the optimal point
! * creates and saves the grids for it
! * only global variables are kept here
! * the most part of the code is universal, and shared by many such modules

module uTMDPDF_OPE
use aTMDe_Numerics
use IntegrationRoutines
use IO_functions
use QCDinput
implicit none

!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v3.00"
character (len=11),parameter :: moduleName="uTMDPDF_OPE"
!Last appropriate version of constants-file
integer,parameter::inputver=0

!--- general
logical:: started=.false.
!! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
integer::outputLevel=2
!! variable that count number of WRNING mesagges. In order not to spam too much
integer::messageTrigger=6

!!! Perturbative order
integer :: orderMain=2 !! LO=0, NLO=1,...

!!! X-Grid parameters
!! over x: i=0...Nx, x_0=xMin
real(dp) :: xMin=0.00001_dp !!! min x 
integer :: Nx=200 !!! number of points in grid
real(dp) :: DeltaX !!! increment of grid
integer :: KminX=0
integer :: KmaxX=5 !!! parameters of range of intepolation

! real(dp) :: xMin=0.1_dp !!! min x 
! integer :: Nx=10 !!! number of points in grid
! real(dp) :: DeltaX !!! increment of grid

!! Numerical parameters
real(dp) :: toleranceINT=1d-8  !!! tolerance for numerical integration

!!! Coefficient lists
integer,parameter::parametrizationLength=37
!! { Log[1-x], log[1-x]^2, log[1-x]^3, log[1-x]^4, log[1-x]^5  !exact
!!   1/x, log[x]/x, Log[x]^2/x  !exact
!! Log[x], log[x]^2, Log[x]^3, log[x]^4, Log[x]^5 !exact
!! T0,...,T23 (Chebyshev polynomials) }
!! The Lmu^2 part is exact the later parts are fitted, but exact if posible (e.g. Lmu and Nf parts for q->q)
!!!!! TO DO: UPDATE TO EXACT VALUES [use 1809.07084]!!


!!--------------------------------------Public interface-----------------------------------------
public::uTMDPDF_OPE_Initialize

!!!!!!----FOR TEST
public::XatNode,invX,NodeForX,Winterpolator!,Tmatrix,TmatrixElement

contains

!! Coefficient function
INCLUDE 'Code/uTMDPDF/coeffFunc-new2.f90'
!INCLUDE 'Code/uTMDPDF/coeffFunc-new.f90'

!! X-grid routines
INCLUDE 'Code/Twist2/Twist2Xgrid.f90'
!! Mellin convolution matrix
!INCLUDE 'Code/Twist2/Twist2MatrixT.f90'


function uTMDPDF_IsInitialized()
    logical::uTMDPDF_IsInitialized
    uTMDPDF_IsInitialized=started
end function uTMDPDF_IsInitialized

!! Initialization of the package
subroutine uTMDPDF_OPE_Initialize(file,prefix)
    character(len=*),intent(in)::file
    character(len=*),intent(in),optional::prefix

    if(started) return
    
    !!!! HERE THE INTIALISATION ROTINE
    
    !!! function to initialze the Xgrid
    call XGrid_Initialize()
    
end subroutine uTMDPDF_OPE_Initialize

end module uTMDPDF_OPE
