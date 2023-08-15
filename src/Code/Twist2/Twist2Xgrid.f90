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
!       Winterpolator(x,i)   - interpolation function such that f(x)=sum_i W(x,i) f(x_i)
!	
!	v.3.00 Created (AV, 18.07.2023)
!
!				A.Vladimirov (18.07.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! sets global variables for the X-griding
subroutine XGrid_Initialize()
    DeltaX=acosh(1_dp/xMin)/Nx
    !DeltaX=-LOG10(xMin)/Nx
end subroutine XGrid_Initialize

!!!! The algorithm is similar to QCD num.
!!!! I.e. it computes the convolution matrices

!!!! Value of X at the node i
pure function XatNode(i)
    real(dp):: XatNode
    integer,intent(in)::i
    
    XatNode=1_dp/cosh(DeltaX*(i-Nx))
    !XatNode=xMin*10**(i*DeltaX)
end function XatNode

!!!! Inverse X
pure function invX(x)
    real(dp):: invX
    real(dp), intent(in):: x
    
    invX=Nx-acosh(1_dp/x)/DeltaX
    !invX=LOG10(x/xMin)/DeltaX
end function invX

!!!! Value of low-grid value for given X
pure function NodeForX(x)
    integer:: NodeForX
    real(dp), intent(in):: x
    
    NodeForX=INT(invX(x))
    
end function NodeForX

!!!! Lagrange polynomial for interpolation
!!!! x- variable,
!!!! j- node to intepolate x->[x_j,x_j+1]
!!!! k=[kMin,...kMax] is the polynomial for y_i+k
!!!! If the segment runs out of grid the corresponding contribution is eliminated
!!!!    it corresponds to a reduced interpolation order.
pure function LagrangeP(x,j,k)
    real(dp):: LagrangeP
    real(dp), intent(in):: x
    integer, intent(in):: j,k
    real(dp)::lx
    integer::l
    
    lx=invX(x)-j
    if(k<KminX .and. k>KmaxX) then
        LagrangeP=0_dp
    else
        LagrangeP=1_dp
        do l=KminX,KmaxX
            !!!  0<= segment <=Nx is a check for segment being in the grid
            if(k/=l .and. j+l<=Nx .and. j+l>=0) then
            LagrangeP=LagrangeP*(lx-l)/(k-l)
            end if
        end do
    end if
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
    if(k<kminX .or. k>kMaxX) then
        Winterpolator=0_dp
    else
        Winterpolator=LagrangeP(x,j,k)
    end if

end function Winterpolator
