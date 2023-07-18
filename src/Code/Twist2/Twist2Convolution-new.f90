!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	This file contains the part of the code, which is common for all TMD-OPE modules 
!       that operates at twist-2. It is inclucded (as a text).
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!	This part is devoted to the calculation of Mellin convolution
!	
!	v.3.00 Created (AV, 18.07.2023)
!
!				A.Vladimirov (18.07.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! The algorithm is similar to QCD num.
!!!! I.e. it computes the convolution matrices

!!!! Value of X at the node i
pure function XatNode(i)
    real(dp):: XatNode
    integer,intent(in)::i
    
    XatNode=10**(-Bx+i*stepX) 
end function XatNode

!!!! Value of low-grid value for given X
pure function NodeForX(x)
    integer:: NodeForX
    real(dp), intent(in):: x
    
    NodeForX=INT((log10(x)+Bx)/stepX)
    
end function NodeForX

!!!! Lagrange polynomial for interpolation
!!!! x- variable,
!!!! j- node to intepolate x->[x_i,x_i+1]
!!!! k=-1,0,1,2 is the polynomial for y_i+k
!!!! The polynomial is cubic, utmost sectors are quadratic
pure function LagrangeP(x,j,k)
    real(dp):: LagrangeP
    real(dp), intent(in):: x
    integer, intent(in):: j,k
    real(dp)::lx
        
    if(j==0) then
        SELECT CASE(k) 
            CASE(0)
                lx=LOG10(x/XatNode(0))/stepX
                LagrangeP=(1-lx)*(1-lx/2)
            CASE(1)
                lx=LOG10(x/XatNode(1))/stepX
                LagrangeP=(1-lx)*(1+lx)
            CASE(2)
                lx=LOG10(x/XatNode(2))/stepX
                LagrangeP=(1+lx)*(1+lx/2)
            CASE DEFAULT
                LagrangeP=0_dp
        END SELECT
    else if(j==Nx-1) then
        SELECT CASE(k) 
            CASE(-1)
                lx=LOG10(x/XatNode(Nx-2))/stepX
                LagrangeP=(1-lx)*(1-lx/2)
            CASE(0)
                lx=LOG10(x/XatNode(Nx-1))/stepX
                LagrangeP=(1-lx)*(1+lx)
            CASE(1)
                lx=LOG10(x/XatNode(Nx))/stepX
                LagrangeP=(1+lx)*(1+lx/2)
            CASE DEFAULT
                LagrangeP=0_dp
        END SELECT
    else
        SELECT CASE(k)
            CASE(-1)
                lx=LOG10(x/XatNode(j-1))/stepX
                LagrangeP=(1-lx)*(1-lx/2)*(1-lx/3)
            CASE(0)
                lx=LOG10(x/XatNode(j))/stepX
                LagrangeP=(1-lx)*(1-lx/2)*(1+lx)
            CASE(1)
                lx=LOG10(x/XatNode(j+1))/stepX
                LagrangeP=(1-lx)*(1+lx/2)*(1+lx)
            CASE(2)
                lx=LOG10(x/XatNode(j+2))/stepX
                LagrangeP=(1+lx)*(1+lx/2)*(1+lx/3)
            CASE DEFAULT
                LagrangeP=0_dp
        END SELECT
    end if
end function LagrangeP
