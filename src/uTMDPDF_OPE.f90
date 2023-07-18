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

INCLUDE 'Tables/G7K15.f90'

!!! Perturbative order
integer :: orderMain=2 !! LO=0, NLO=1,...

!!! Grid parameters
!! over x: x_i=10^(-Bx+Delta*i), i=0...Nx
real(dp) :: Bx=5._dp !!! min x is 10^{-Bx}
integer :: Nx=200 !!! number of points in grid
real(dp) :: stepX=0.025_dp !!! increment of grid

!!!!!!----FOR TEST
public::XatNode,NodeForX,LagrangeP

contains

INCLUDE 'Code/Twist2/Twist2Convolution-new.f90'

end module uTMDPDF_OPE
