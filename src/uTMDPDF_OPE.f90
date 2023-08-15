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
use uTMDPDF_model_new
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
!! variable that count number of WARNING mesagges. In order not to spam too much
integer::messageTrigger=6

!!!------------------------- PARAMETERS DEFINED IN THE INI-file--------------

!!! Perturbative order
integer :: orderMain=3 !! LO=0, NLO=1,...

!!!! X-Grid parameters
!! over x: i=0...Nx, x_0=xMin
real(dp) :: xMin=0.00001_dp !!! min x 
integer :: Nx=200 !!! number of points in grid
real(dp) :: DeltaX !!! increment of grid
integer :: KminX=-1
integer :: KmaxX=2 !!! parameters of range of intepolation

!!!! Numerical parameters
real(dp) :: toleranceINT=1d-6  !!! tolerance for numerical integration
real(dp) :: bFREEZE=1d-6       !!! small value of b at which b freesed

logical(dp) :: IsMuYdependent  !!! if mu is y independent, computation is much(!) faster

!!!------------------------- HARD-CODED PARAMETERS ----------------------
!!! Coefficient lists
integer,parameter::parametrizationLength=37
!! { Log[1-x], log[1-x]^2, log[1-x]^3, log[1-x]^4, log[1-x]^5  !exact
!!   1/x, log[x]/x, Log[x]^2/x  !exact
!! Log[x], log[x]^2, Log[x]^3, log[x]^4, Log[x]^5 !exact
!! T0,...,T23 (Chebyshev polynomials) }
!! The Lmu^2 part is exact the later parts are fitted, but exact if posible (e.g. Lmu and Nf parts for q->q)
!!!!! TO DO: UPDATE TO EXACT VALUES [use 1809.07084]!!

!!!------------------------- DYNAMICAL-GLOBAL PARAMETERS -------------------
real(dp) :: c4_global=1_dp  !!! scale variation parameter

!!--------------------------------------Public interface-----------------------------------------
public::uTMDPDF_OPE_Initialize

!!!!!!----FOR TEST
public::XatNode,invX,NodeForX,Winterpolator,CxF_compute!,Tmatrix,TmatrixElement

contains

!! Coefficient function
INCLUDE 'Code/uTMDPDF/coeffFunc-new2.f90'
!INCLUDE 'Code/uTMDPDF/coeffFunc-new.f90'

!! X-grid routines
INCLUDE 'Code/Twist2/Twist2Xgrid.f90'
!! Mellin convolution matrix
INCLUDE 'Code/Twist2/Twist2Convolution-new.f90'
!INCLUDE 'Code/Twist2/Twist2MatrixT.f90'


function uTMDPDF_IsInitialized()
    logical::uTMDPDF_IsInitialized
    uTMDPDF_IsInitialized=started
end function uTMDPDF_IsInitialized

!! Initialization of the package
subroutine uTMDPDF_OPE_Initialize(file,prefix)
    character(len=*),intent(in)::file
    character(len=*),intent(in),optional::prefix
    integer::i

    if(started) return
    
    !!!! HERE THE INTIALISATION ROTINE
    
    !!! function to initialze the Xgrid
    call XGrid_Initialize()
    
    !!! _OPE_model initialisation
    !call ModelInitialization()
    
    !!!!!!!Checking the x-dependance of muOPE
    IsMuYdependent=testMU()

    if(IsMuYdependent) then
        if(outputLevel>2) write(*,*) trim(moduleName)//': mu OPE is dependent on x'
    else
        if(outputLevel>2) write(*,*) trim(moduleName)//': mu OPE is independent on x'
    end if
    
end subroutine uTMDPDF_OPE_Initialize

!!!!array of x times PDF(x,Q) for hadron 'hadron'
!!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function xf(x,Q,hadron)
    real(dp) :: x,Q
    integer:: hadron
    real(dp), dimension(-5:5):: xf
    
    xf=xPDF(x,Q,hadron)
    
end function xf

end module uTMDPDF_OPE
