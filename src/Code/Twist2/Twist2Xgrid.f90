!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	This file contains the part of the code, which is common for all TMD-OPE modules 
!       that operates at twist-2. It is inclucded (as a text).
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!       The grid is hard coded to be cubic forward KK=[0,1,2,3]
!	This part is devoted to the common functions for X-grid
!       XGrid_Initialize()
!       real XatNode(i)      - gives x_i
!       real invX(x)         - gives i such that x_i=x
!       integer NodeForX(x)  - gives i such that x_i<x<x_{i+1}
!       Winterpolator(x,i)   - interpolation function such that f(x)=sum_i W(x,i) f(x_i)
!	
!	v.3.00 Created (AV, 18.07.2023)
!
!				A.Vladimirov (18.07.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! sets global variables for the X-griding
subroutine XGrid_Initialize()
    DeltaX=acosh(1/xMin)/Nx
end subroutine XGrid_Initialize

!!!! The algorithm is similar to QCD num.
!!!! I.e. it computes the convolution matrices

!!!! Value of X at the node i
pure function XatNode(i)
    real(dp):: XatNode
    integer,intent(in)::i
    
    XatNode=1_dp/cosh(DeltaX*(i-Nx))
end function XatNode

!!!! Inverse X
pure function invX(x)
    real(dp):: invX
    real(dp), intent(in):: x
    
    invX=Nx-acosh(1/x)/DeltaX
    
end function invX

!!!! Value of low-grid value for given X
pure function NodeForX(x)
    integer:: NodeForX
    real(dp), intent(in):: x
    
    NodeForX=INT(invX(x))
    
end function NodeForX

!!!! Lagrange polynomial for interpolation
!!!! x- variable,
!!!! j- node to intepolate x->[x_i,x_i+1]
!!!! k=0,1,2,3 is the polynomial for y_i+k
!!!! The polynomial is cubic forward 
!!!! WARNING!!! there is no check for x and j being in grid
!!!! WARNING!!! cubic forwardness is hard coded. It is used in the constraction of Tmatrix
pure function LagrangeP(x,j,k)
    real(dp):: LagrangeP
    real(dp), intent(in):: x
    integer, intent(in):: j,k
    real(dp)::lx
    
    lx=invX(x)-j    
    SELECT CASE(k)
        CASE(0)
            LagrangeP=(1-lx)*(2-lx)*(3-lx)/6
        CASE(1)
            LagrangeP=lx*(2-lx)*(3-lx)/2
        CASE(2)
            LagrangeP=-lx*(1-lx)*(3-lx)/2
        CASE(3)
            LagrangeP=lx*(1-lx)*(2-lx)/6
        CASE DEFAULT
            LagrangeP=0_dp
    END SELECT
end function LagrangeP


!!! The intepolation function, which gives x->i
!!! i.e. A PDF is presented as f(x)=W_i(x)f(x_i)
!!! The function is W=O(x_{i-3}<x<x_{i-2})P(x,i-3,3)+...+O(x_{i}<x<x_{i+1})P(x,i,0)
!!! x, and i are natural arguments.
!!! If x outside of grid ->0
pure function Winterpolator(x,i)
    real(dp)::Winterpolator
    real(dp),intent(in)::x
    integer,intent(in)::i
    integer::j,k
    
    !!! first get the node for x
    j=NodeForX(x)
    !!! if node out side of the grid. >>0
    if(j<0 .or. j>Nx-1) then
        Winterpolator=0_dp
        return
    end if
    !!! the difference is the shift of interpolator
    k=i-j
    Winterpolator=LagrangeP(x,j,k)

end function Winterpolator
