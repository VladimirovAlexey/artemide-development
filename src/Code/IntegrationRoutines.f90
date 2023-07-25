!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.04
!
!	The module that contains standard integration routines (mainly used for bin integration)
!       the main difference from library routines is the absence of checks (boundaries, etc.), 
!           the checks are done in the artemide main code.
!
!   Currently it incorporates:
!       S5: 5-point Simpsons
!       SN: N-point Simpsons (N=even)
!       SA: Adaptive Simpsons
!       G7: 7-point Gauss
!       K15: 15-point Kronrod
!       GK: Adaptive Gauss Kronrod 7/15
!
!				A.Vladimirov (17.04.2020)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module IntegrationRoutines
use aTMDe_numerics
implicit none

private
INCLUDE '../Tables/G7K15.f90'

public::Integrate_S5,Integrate_SN,Integrate_SA
public::Integrate_G7,Integrate_K15,Integrate_GK
public::Integrate_GK_array5

!!! this is interface for function (-5:5) in the integration
abstract interface 
    function func_array5(x)
        import::dp
        real(dp),dimension(-5:5) :: func_array5
        real(dp), intent(in) ::x
    end function func_array5
end interface
  
contains


!!! Simpson by 5 points
!!! Use it for estimations only!
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
function Integrate_S5(f,xMin,xMax)
    real(dp)::f,Integrate_S5
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,x2,x3,x4,f1,f2,f3,f4,f5

    delta=xMax-xMin
    x2=xMin+delta/4._dp
    x3=xMin+delta/2._dp
    x4=xMax-delta/4._dp

    f1=f(xMin)
    f2=f(x2)
    f3=f(x3)
    f4=f(x4)
    f5=f(xMax)
        
    Integrate_S5=delta*(f1+4._dp*f2+2._dp*f3+4._dp*f4+f5)/12._dp
    
end function Integrate_S5


!!! Simpson by N points
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
!!! N is an even integer N>2
function Integrate_SN(f,xMin,xMax,N)
    real(dp)::f,Integrate_SN
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,inter,xCur
    integer::N,i

    !!!!!!!!!!!!!!!!!!!fixed number Simpsons
    delta=(xMax-xMin)/N
    inter=f(xMin)   !!!! first term
    
    !!!! even terms
    do i=1,N-1,2
        xCur=xMin+i*delta
        inter=inter+4._dp*f(xCur)    
    end do
    !!!! odd term
    do i=2,N-2,2
        xCur=xMin+i*delta
        inter=inter+2._dp*f(xCur)    
    end do
    
    inter=inter+f(xMax)!!!! last term
        
    Integrate_SN=delta/3._dp*inter
    
end function Integrate_SN

!!! Simpson adaptive
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
!!! tolerance is relative (it is wieghted by approximate value of integral)
function Integrate_SA(f,xMin,xMax,tolerance)
    real(dp)::f,Integrate_SA,tolerance,eps
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,x2,x3,x4,f1,f2,f3,f4,f5

    
    delta=xMax-xMin
    x2=xMin+delta/4._dp
    x3=xMin+delta/2._dp
    x4=xMax-delta/4._dp

    f1=f(xMin)
    f2=f(x2)
    f3=f(x3)
    f4=f(x4)
    f5=f(xMax)
    
    !!! the error parameter is weighted with the approximate integral size
    eps=tolerance*abs(delta*(f1+4._dp*f2+2._dp*f3+4._dp*f4+f5)/12._dp)
        
    Integrate_SA=SA_Rec(f,xMin,x2,x3,f1,f2,f3,eps)+SA_Rec(f,x3,x4,xMax,f3,f4,f5,eps)
    
end function Integrate_SA

recursive function SA_Rec(f,x1,x3,x5,f1,f3,f5,eps) result(res)
    real(dp)::f,x1,x2,x3,x4,x5
    real(dp)::f1,f2,f3,f4,f5,eps,res
    real(dp)::value15,value135
    
    x2=(x1+x3)/2._dp
    f2=f(x2)
    x4=(x3+x5)/2._dp
    f4=f(x4)
    
    value15=(x5-x1)*(f1+4._dp*f3+f5)/6._dp
    value135=(x5-x1)*(f1+4._dp*f2+2._dp*f3+4._dp*f4+f5)/12._dp
        
    If(ABS(value135-value15)>eps) then
        res=SA_Rec(f,x1,x2,x3,f1,f2,f3,eps)+SA_Rec(f,x3,x4,x5,f3,f4,f5,eps)
    else
        res=value135
    end if
    
end function SA_Rec

!!! Gauss 7-points
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
function Integrate_G7(f,xMin,xMax)
    real(dp)::f,Integrate_G7
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,inter
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp
    
    inter=0._dp
    do i=1,7
        inter=inter+Wi_g77(i)*f(Xi_g7(i)*delta+av)
    end do
        
    Integrate_G7=delta*inter
    
end function Integrate_G7

!!! Kronrod 15-points
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
function Integrate_K15(f,xMin,xMax)
    real(dp)::f,Integrate_K15
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,inter
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp
    
    inter=0._dp
    do i=1,15
        inter=inter+Wi_k15(i)*f(Xi_k15(i)*delta+av)
    end do
        
    Integrate_K15=delta*inter
    
end function Integrate_K15


!!! Gauss-Kronrod 7/15 adaptive
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
!!! tolerance is relative (it is wieghted by approximate value of integral)
function Integrate_GK(f,xMin,xMax,tolerance)
    real(dp)::f,Integrate_GK
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,g7,k15,eps,tolerance,fI
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp
    
    g7=0._dp
    k15=0._dp
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

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp
    
    g7=0._dp
    k15=0._dp
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

!!! Gauss-Kronrod 7/15 adaptive for array
!!! f::  array(-5:5) function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
!!! tolerance is relative (it is wieghted by approximate value of integral)
!!! the integral is considered convergent once all components are convergent
function Integrate_GK_array5(f,xMin,xMax,tolerance)
    real(dp),dimension(-5:5)::Integrate_GK_array5
    procedure(func_array5)::f
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,tolerance
    real(dp), dimension(-5:5)::fI,g7,k15,eps
    integer::i
    logical::ISconvergent
    
    
    
    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp
    g7=0._dp
    k15=0._dp
    do i=1,15
        fI=f(Xi_k15(i)*delta+av)
        g7=g7+Wi_g7(i)*fI
        k15=k15+Wi_k15(i)*fI
    end do
    
    !!! check convergence
    eps=delta*abs(k15)*tolerance  
    
    ISconvergent=.true.
    do i=-5,5
       if(delta*abs(k15(i)-g7(i))>eps(i)) then
        ISconvergent=.false.
        exit
       end if
    end do
    
    if(ISconvergent) then
        Integrate_GK_array5=delta*k15                
    else
        Integrate_GK_array5=GK_array5_Rec(f,xMin,av,eps)+GK_array5_Rec(f,av,xMax,eps)
    end if
    
end function Integrate_GK_array5

recursive function GK_array5_Rec(f,xMin,xMax,eps) result(res)
    real(dp),dimension(-5:5)::res
    procedure(func_array5)::f
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av
    real(dp),dimension(-5:5)::g7,k15,eps,fI
    integer::i
    logical::ISconvergent
    
    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp

    g7=0._dp
    k15=0._dp
    do i=1,15
        fI=f(Xi_k15(i)*delta+av)
        g7=g7+Wi_g7(i)*fI
        k15=k15+Wi_k15(i)*fI
    end do
    
    !!! check convergence
    ISconvergent=.true.
    do i=-5,5
       if(delta*abs(k15(i)-g7(i))>eps(i)) then
        ISconvergent=.false.
        exit
       end if
    end do
    
    if(ISconvergent) then
        res=delta*k15                
    else
        res=GK_array5_Rec(f,xMin,av,eps)+GK_array5_Rec(f,av,xMax,eps)
    end if
    
end function GK_array5_Rec

end module IntegrationRoutines
