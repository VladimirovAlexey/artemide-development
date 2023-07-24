!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model parameters for unpolarized TMDPDF  small-b OPE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module uTMDPDF_OPE_model
use aTMDe_Numerics
use IO_functions
implicit none

private

!!!!!------------------------------------------------------------------------------------
!!!!! These functions MUST defined in module  !!
!!!!! 0) Subroutine which is called at the start. No arguments
public:: ModelInitialization
!!!!! 1) Function which returns the value of b used as argument of convolution integrals
!!!!!    arg=(bT,x,y) with br=transverse distance(real_dp), x being Bjorken x (real_dp), and y being convolution variable (real_dp)
real(dp),public:: bSTAR
!!!!! 2) Function which returns the scale of matching (OPE scale)
!!!!!    arg=(bT,x,y,c4) with bT, x, y same as in bSTAR, and c4(real_dp) is scale variation parameter
real(dp),public:: muOPE
!!!!!------------------------------------------------------------------------------------


contains  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER DEFINED FUNCTIONS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!! Write nessecery model intitialization.
subroutine ModelInitialization()    
    write(*,*) color(">>>  The model for uTMDPDF OPE for ART23   <<<",c_cyan)    
end subroutine ModelInitialization

!!!! This is the function b* that enters the logarithms of coefficient function
!!!! at small-b it should be ~b to match the collinear regime
!!!! at large-b it is a part of model
!!!! x -- is the global x for TMDPDF,
!!!! y -- is the convolution variable in the definition \int dy/y C(y) PDF(x/y)
pure function bSTAR(bT,x,y)
    real(dp),intent(in)::bT,x,y

    bSTAR=bT/sqrt(1d0+(bT/500d0)**2)
end function bSTAR
  
!!!!This function is the mu(x,b), which is used inside the OPE
!!!! x -- is the global x for TMDPDF,
!!!! y -- is the convolution variable in the definition \int dy/y C(y) PDF(x/y)
!!!! c4-- is the scale variation variable
pure function muOPE(bt,x,y,c4)
    real(dp),intent(in)::bt,x,y,c4

    muOPE=C0_const*c4/bT+2d0
    
    if(muOPE>1000d0) then
        muOPE=1000d0
    end if
end function muOPE
  
end module uTMDPDF_OPE_model
