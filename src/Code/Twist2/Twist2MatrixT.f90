!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	This file contains the part of the code, which is common for all TMD-OPE modules 
!       that operates at twist-2. It is inclucded (as a text).
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!	This part is devoted to the calculation of Mellin convolution between
!       interpolation function and the coefficient function.
!       Tmatrix(flav,n,k,Nf) = real(0:Nx,0:Nx)
!       int flav: flavor combination 1=qq,2=qg,3=gq,4=gg,5=qqb,6=qqp
!       int n, k: as^n L^k
!       real Nf: number of flavors
!
!       Tmatrix2 = Tmatrix but it uses symmetries and MUCH faster
!	
!	v.3.00 Created (AV, 18.07.2023)
!
!				A.Vladimirov (18.07.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!! The matrix T is defined as Mellin convolution of [C W_j]_i
!!!! The parameter flav is unversal flavor parameter
!!!! i,j is the element of the matrix
!!!! n,k, Nf are standard variables
function TmatrixElement(flav,i,j,n,k,Nf)
    real(dp)::TmatrixElement
    integer, intent(in)::i,j,n,k,flav
    real(dp),intent(in)::Nf
    real(dp)::inter,vMin,vMax,xCurrent,xi
    real(dp),dimension(1:3)::singC
    
    inter=0_dp
    
    xi=XatNode(i)  
    
    if(i==j) then             
        if(flav==1) then
            !!! delta-function contribution
            inter=C_q_q_delta(n,k,Nf)            
            !!! The remnant [0,x_i] from the (..)_+ integration
            singC=C_q_q_plus(n,k,Nf)
            xCurrent=xi/XatNode(j+1)
            inter=inter+Log(1-xCurrent)*singC(1)
            inter=inter+Log(1-xCurrent)**2*singC(2)/2_dp
            inter=inter+Log(1-xCurrent)**3*singC(3)/3_dp
            
        else if(flav==4) then
            inter=C_g_g_delta(n,k,Nf)
            
            !!! The remnant [0,x_i] from the (..)_+ integration
            singC=C_g_g_plus(n,k,Nf)
            xCurrent=xi/XatNode(j+1)
            inter=inter+Log(1-xCurrent)*singC(1)
            inter=inter+Log(1-xCurrent)**2*singC(2)/2_dp
            inter=inter+Log(1-xCurrent)**3*singC(3)/3_dp
        end if
    end if        
    
    if(i<j+1) then
        !! compute boundaries
        if(j==Nx) then
            vMin=xi
        else
            vMin=xi/XatNode(j+1)
        end if
        
        if(i<j-3) then
            vMax=xi/XatNode(j-3)
        else
            vMax=1_dp
        end if
        
        !! regular part all flavors
        inter=inter+Integrate_GK(FFreg,vMin,vMax,toleranceINT)  
        
        !! (..)_+ part
        if(flav==1 .or. flav==4) then
            !! (..)_+ part only for qq and gg
            if(i==j) then
                !! note that in this case singC is compute earlier
                inter=inter+Integrate_GK(FFplus,vMin,vMax,toleranceINT)
            else
                inter=inter+Integrate_GK(FFplus2,vMin,vMax,toleranceINT)
            end if        
        end if
    end if            
    
    TmatrixElement=inter

    !!!!! these are wrapper functions to pass ther integrand to GK routine
    !!!!! FORTRAN gets only function of one variable
    !!!!! using the trick with internal function one can by-pass this limitation
    contains
    
    !!! regular integrand
    function FFreg(x)
        real(dp)::FFreg
        real(dp),intent(in)::x
        FFreg=CoefREG(flav,x,n,k,Nf)*Winterpolator(xi/x,j)
    end function FFreg
    
    !!! (..)_+ integrand
    function FFplus(x)
        real(dp)::FFplus
        real(dp),intent(in)::x
        FFplus=(singC(1)+singC(2)*Log(1-x)+singC(2)*Log(1-x)**2)/(1-x)*&
            (Winterpolator(xi/x,j)-1)
    end function FFplus
    
    !!! (..)_+ integrand, that does not hit the delta-function
    function FFplus2(x)
        real(dp)::FFplus2
        real(dp),intent(in)::x
        FFplus2=(singC(1)+singC(2)*Log(1-x)+singC(2)*Log(1-x)**2)/(1-x)*&
            Winterpolator(xi/x,j)
    end function FFplus2
    
end function TmatrixElement

!!!! The matrix T is defined as Mellin convolution of [C W_j]_i
!!!! The parameter flav is unversal flavor parameter
function Tmatrix(flav,n,k,Nf)
    real(dp),dimension(0:Nx,0:Nx)::Tmatrix
    integer, intent(in)::n,k,flav
    real(dp),intent(in)::Nf
    integer::i,j
    
    Tmatrix=0_dp
    
    do i=0,Nx-1
    do j=0,Nx        
        Tmatrix(i,j)=TmatrixElement(flav,i,j,n,k,Nf)
    end do
    end do
end function Tmatrix

!!!! The matrix T is defined as Mellin convolution of [C W_j]_i
!!!! The parameter flav is unversal flavor parameter
function Tmatrix2(flav,n,k,Nf)
    real(dp),dimension(0:Nx,0:Nx)::Tmatrix2
    integer, intent(in)::n,k,flav
    real(dp),intent(in)::Nf
    integer::i,j
    real(dp),dimension(0:Nx-3)::matrixSequence
    
    Tmatrix2=0_dp
    
    !!! The first line is unique
    do j=0,Nx        
        Tmatrix2(0,j)=TmatrixElement(flav,0,j,n,k,Nf)
    end do
    
    !!! the 0:n-3 elements of matrix are repeated.
    !!! save them to a separate value
    do j=0,Nx-3        
        matrixSequence(j)=TmatrixElement(flav,1,j,n,k,Nf)
    end do
    
    !!!! 0:N-3 from sequence rest is computed
    do i=1,Nx-1
        Tmatrix2(i,i-1:Nx-3)=matrixSequence(0:Nx-i-2)
        Tmatrix2(i,Nx-2)=TmatrixElement(flav,i,Nx-2,n,k,Nf)
        Tmatrix2(i,Nx-1)=TmatrixElement(flav,i,Nx-1,n,k,Nf)
        Tmatrix2(i,Nx)=TmatrixElement(flav,i,Nx,n,k,Nf)
    end do
end function Tmatrix2
