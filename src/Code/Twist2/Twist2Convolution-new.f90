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

!!!! Lagrange polynomial for interpolation over Log10-grid
!!!! x- variable,
!!!! j- node to intepolate x->[x_i,x_i+1]
!!!! k=-1,0,1,2 is the polynomial for y_i+k
!!!! The polynomial is cubic, utmost sectors are quadratic
!!!! WARNING!!! there is no check for x and j being in grid
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


!!! The intepolation function, which gives x->i
!!! i.e. A PDF is presented as f(x)=W_i(x)f(x_i)
!!! The function is W=O(x_{i-2}<x<x_{i-1})P(x,i-2,2)+...+O(x_{i+1}<x<x_{i+1})P(x,i+1,-1)
!!! x, and i are natural arguments.
!!! If x outside of grid ->0
pure function Winterpolator(x,i)
    real(dp)::Winterpolator
    real(dp),intent(in)::x
    integer,intent(in)::i
    integer::j,k
    
    !!! first get the node for x
    j=NodeForX(x)
    !!! if
    if(j<0 .or. j>Nx-1) then
        Winterpolator=0_dp
        return
    end if
    !!! the difference is the shift of interpolator
    k=i-j
    if(k<-1 .or. k>2) then
        Winterpolator=0_dp
    else
        Winterpolator=LagrangeP(x,j,k)
    end if

end function Winterpolator


function Tmatrix(flav,n,k,Nf)
    real(dp),dimension(0:Nx,0:Nx)::Tmatrix
    integer, intent(in)::n,k,flav
    real(dp),intent(in)::Nf
    integer::i,j
    real(dp)::inter,vMin,vMax,xCurrent
    real(dp),dimension(1:3)::singC
    
    Tmatrix=0_dp
    
    do i=0,Nx-1
    do j=0,Nx
        inter=0_dp
        
        if(i==j) then 
            !!! delta-function contribution
            if(flav==1) then
                inter=C_q_q_delta(n,k,Nf)
                
                !!! The remiant [0,x_i] from the (..)_+ integration
                singC=C_q_q_plus(n,k,Nf)
                xCurrent=XatNode(i)
                inter=inter+Log(1-xCurrent)*singC(1)
                inter=inter+Log(1-xCurrent)**2*singC(2)/2_dp
                inter=inter+Log(1-xCurrent)**3*singC(3)/3_dp
            else if(flav==4) then
                inter=C_g_g_delta(n,k,Nf)
                
                !!! The remiant [0,x_i] from the (..)_+ integration
                singC=C_g_g_plus(n,k,Nf)
                xCurrent=XatNode(i)
                inter=inter+Log(1-xCurrent)*singC(1)
                inter=inter+Log(1-xCurrent)**2*singC(2)/2_dp
                inter=inter+Log(1-xCurrent)**3*singC(3)/3_dp
            end if
            
            
        end if
        
        !! regular part
        if(i<j+2) then
            vMin=XatNode(i)/XatNode(j+2)
            if(i<j-2) then
                vMax=XatNode(i)/XatNode(j-2)
            else
                vMax=1_dp
            end if
            inter=inter+integrateREG_q_q(flav,i,j,n,k,Nf,vMin,vMax)
            if(flav==1 .or. flav==4) then
            inter=inter+integrateSING_q_q(flav,i,j,n,k,Nf,vMin,vMax)
            end if
        end if            
    
    Tmatrix(i,j)=inter
    end do
    end do
end function Tmatrix

!!!!! this is a wrapper functions to pass ther integrand to GK routine
function integrateREG_q_q(flav,i,j,n,k,Nf,vMin,vMax)
    integer, intent(in)::i,j,n,k,flav
    real(dp),intent(in)::Nf,vMin,vMax
    real(dp)::integrateREG_q_q
    real(dp)::interX
    
    interX=XatNode(i)  
    
    integrateREG_q_q=Integrate_GK(FF,vMin,vMax,toleranceINT)
    
    contains
    
    function FF(x)
        real(dp)::FF
        real(dp),intent(in)::x
        FF=CoefREG(flav,x,n,k,Nf)*Winterpolator(interX/x,j)
    end function FF
end function integrateREG_q_q

!!!!! this is a wrapper function to pass ther integrand to GK routine
function integrateSING_q_q(flav,i,j,n,k,Nf,vMin,vMax)
    integer, intent(in)::i,j,n,k,flav
    real(dp),intent(in)::Nf,vMin,vMax
    real(dp)::integrateSING_q_q
    real(dp)::interX
    real(dp),dimension(1:3)::singC
    
    interX=XatNode(i)      
    singC=CoefSING(flav,n,k,Nf)
    
    if(i==j) then
    integrateSING_q_q=Integrate_GK(FF,vMin,vMax,toleranceINT)
    else
    integrateSING_q_q=Integrate_GK(FF2,vMin,vMax,toleranceINT)
    end if
    
    contains
    
    function FF(x)
        real(dp)::FF
        real(dp),intent(in)::x
        FF=(singC(1)+singC(2)*Log(1-x)+singC(2)*Log(1-x)**2)/(1-x)*&
            (Winterpolator(interX/x,j)-1)
    end function FF
    
    function FF2(x)
        real(dp)::FF2
        real(dp),intent(in)::x
        FF2=(singC(1)+singC(2)*Log(1-x)+singC(2)*Log(1-x)**2)/(1-x)*&
            Winterpolator(interX/x,j)
    end function FF2
    
end function integrateSING_q_q

!!! ----------------------------------------------------------------
!!! Gauss-Kronrod 7/15 adaptive
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
!!! tolerance is relative (it is wieghted by approximate value of integral)
function Integrate_GK(f,xMin,xMax,tolerance)
    real(dp)::f,Integrate_GK
    real(dp),intent(in)::xMin,xMax,tolerance
    real(dp)::delta,av,g7,k15,eps,fI
    integer::i

    delta=(xMax-xMin)/2_dp
    av=(xMax+xMin)/2_dp
    
    g7=0_dp
    k15=0_dp
    do i=1,15
        fI=f(Xi_k15(i)*delta+av)
        g7=g7+Wi_g7(i)*fI
        k15=k15+Wi_k15(i)*fI
    end do
    
    eps=delta*abs(k15)*tolerance
    
    if(delta*abs(k15-g7)>eps) then
        Integrate_GK=GK_Rec(f,xMin,av,eps)+GK_Rec(f,av,xMax,eps)
    else
        Integrate_GK=delta*k15
    end if
    
end function Integrate_GK

recursive function GK_Rec(f,xMin,xMax,eps) result(res)
    real(dp)::f,res
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,g7,k15,eps,fI
    integer::i

    delta=(xMax-xMin)/2_dp
    av=(xMax+xMin)/2_dp
    
    g7=0_dp
    k15=0_dp
    do i=1,15
        fI=f(Xi_k15(i)*delta+av)
        g7=g7+Wi_g7(i)*fI
        k15=k15+Wi_k15(i)*fI
    end do
    
    if(delta*abs(k15-g7)>eps) then
        res=GK_Rec(f,xMin,av,eps)+GK_Rec(f,av,xMax,eps)
    else
        res=delta*k15
    end if
    
end function GK_Rec

