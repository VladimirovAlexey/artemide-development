!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	This file contains the part of the code, which is common for all TMD-OPE modules 
!       that operates at twist-2. It is inclucded (as a text).
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!       The interpolation type is defined by K, 
!       as the segments runs out of grid, the interpolation order is reduced 
!	This part is devoted to the common functions for X-grid
!       XGrid_Initialize()
!       real XatNode(i)      - gives x_i
!       real invX(x)         - gives i such that x_i=x
!       integer NodeForX(x)  - gives i such that x_i<x<x_{i+1}
!	
!	v.3.00 Created (AV, 18.07.2023)
!
!				A.Vladimirov (18.07.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! sets global variables for the X-griding
subroutine XGrid_Initialize()
    DeltaX=(acosh(1._dp/xMin)**(1._dp/parX))/Nx
end subroutine XGrid_Initialize

!!!! Value of X at the node i
pure function XatNode(i)
    real(dp):: XatNode
    integer,intent(in)::i
    
    XatNode=1._dp/cosh((DeltaX*(Nx-i))**parX)
end function XatNode

!!!! Value of X at the node i(real)
pure function XatNode_real(i)
    real(dp):: XatNode_real
    real(dp),intent(in)::i

    XatNode_real=1._dp/cosh((DeltaX*(Nx-i))**parX)
end function XatNode_real

!!!! Inverse X
pure function invX(x)
    real(dp):: invX
    real(dp), intent(in):: x
    
    invX=Nx-(acosh(1._dp/x)**(1._dp/parX))/DeltaX
end function invX

!!!! Value of low-grid value for given X
pure function NodeForX(x)
    integer:: NodeForX
    real(dp), intent(in):: x
    
    NodeForX=INT(invX(x))
end function NodeForX
