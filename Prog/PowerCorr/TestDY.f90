module ONE
use aTMDe_Numerics
use aTMDe_control
use EWinput
use QCDinput
use TMDs_inKT
implicit none

private

INCLUDE '../../src/Tables/G7K15.f90'

integer::M1

public:: SetM1,integralKK_GK

contains

subroutine SetM1(m)
    integer::m
    M1=m
end subroutine SetM1

! function uTMDPDF_kT_5(x,kT,a,b,c)
!     real*8::x,kT,a,b,lll
!     integer::c
!     real*8::uTMDPDF_kT_5(-5:5)
!     lll=(1d0+0.2d0*sqrt(x))*(1d0-x)**(0.3d0)/x*(exp(-0.2*kT**2)+2d0/(1d0+kT**2))
!     
!     uTMDPDF_kT_5=lll*(/1d0,1d0,1d0,1d0,1d0,0d0,1d0,1d0,1d0,1d0,1d0/)
! 
! end function uTMDPDF_kT_5

!!! expressions for xi_{1,2}
function xi1(K,dK,tt2)
    real*8::xi1,K,dK,tt2
    xi1=(tt2+K*dK+sqrt((tt2-K**2)*(tt2-dK**2)))/2d0/tt2        
    
    if(ISNAN(xi1)) then
        write(*,*) tt2,K**2,dK**2
    end if
end function xi1

function xi2(K,dK,tt2)
    real*8::xi2,K,dK,tt2
    xi2=(tt2-K*dK+sqrt((tt2-K**2)*(tt2-dK**2)))/2d0/tt2
end function xi2

!!!! This in the integrand without singular roots
function FF(x1,x2,Q2,tau2,qT,K,dK)
    real*8::FF,x1,x2,Q2,K,dK,tau2,qT
    real*8::KK,dKK,xx1,xx2
    real*8::FA(-5:5),FB(-5:5),FAB(-5:5),TT,TT0
    
    KK=K**2
    dKK=dK**2
    
    if(M1==1 .or. M1==3) then
        xx1=xi1(K,dK,tau2)
        xx2=xi2(K,dK,tau2)        
        
        FA=uTMDPDF_kT_5(x1*xx1,(K+dK)/2d0,sqrt(Q2),Q2,1)
        FB=uTMDPDF_kT_5(x2*xx2,(K-dK)/2d0,sqrt(Q2),Q2,1)
        
        FAB=FA*FB
	
        TT=FAB(1)/9.d0&
            +FAB(2)*4.d0/9.d0&
            +FAB(3)/9.d0&
            +FAB(4)*4d0/9.d0&
            +FAB(5)/9d0&
            +FAB(-1)/9.d0&
            +FAB(-2)*4.d0/9.d0&
            +FAB(-3)/9.d0&
            +FAB(-4)*4d0/9.d0&
            +FAB(-5)/9d0
        
        TT=(2d0*tau2-KK-dKK+2d0*sqrt((tau2-KK)*(tau2-dKK)))/4d0/(sqrt(tau2-dKK)*sqrt(tau2-KK))/(xx1*xx2*(1+qT**2/Q2))*TT
	 else
        TT=0d0	 
	 end if
	 
	if(M1==1 .or. M1==2 .or. M1==4) then
        FA=uTMDPDF_kT_5(x1,(K+dK)/2d0,sqrt(Q2),Q2,1)
        FB=uTMDPDF_kT_5(x2,(K-dK)/2d0,sqrt(Q2),Q2,1)    
        
        FAB=FA*FB
	
        TT0=FAB(1)/9.d0&
            +FAB(2)*4.d0/9.d0&
            +FAB(3)/9.d0&
            +FAB(4)*4d0/9.d0&
            +FAB(5)/9d0&
            +FAB(-1)/9.d0&
            +FAB(-2)*4.d0/9.d0&
            +FAB(-3)/9.d0&
            +FAB(-4)*4d0/9.d0&
            +FAB(-5)/9d0
        
    end if
  
    if(M1==1) then
        FF=TT-TT0
    else if(M1==2) then
        FF=TT0
    else if(M1==3) then
        FF=TT
    end if    
end function FF

!!! to integrate over dK I split the integral {-qT,0,qT} and in each sector subtract the singularity at qT
!!! the formula is 
!!! \int_{-qT}^{qT} f(k)/sqrt(qT^2-k^2)dk = \int_0^qT (f(k)/sqrt(qT+k)-f(qT)/sqrt(2 qT))/sqrt(qT-k)+sqrt(2)f(qT)
!!!                                         +\int_{-qT}^0 (f(k)/sqrt(qT-k)-f(-qT)/sqrt(2 qT))/sqrt(qT+k)+sqrt(2)f(-qT)
function integraldK(x1,x2,Q2,tau2,qT,K)
    real*8::integraldK,x1,x2,Q2,tau2,K,qT
    real*8::deltaK,F0,valueMax
    real*8::K2,K3,K4,K5,K6,K7,K8
    real*8::XX1,XX2,XX3,XX4,XX5,XX6,XX7,XX8,XX9
    
!     do i=-10,10
!     write(*,*) "{",qT*i/10d0,",",FF(x1,x2,Q2,tau2,K,qT*i/10d0),"},"
!     end do
    
    !!! k>0 integral
    deltaK=qT
    K2=deltaK/8d0
    K3=deltaK/4d0
    K4=3d0*deltaK/8d0
    K5=deltaK/2d0
    K6=5d0*deltaK/8d0
    K7=3d0*deltaK/4d0
    K8=7d0*deltaK/8d0
    
    
    F0=FF(x1,x2,Q2,tau2,qT,K,qT)/sqrt(2d0*qT)
    
    XX1=(FF(x1,x2,Q2,tau2,qT,K,0d0)/sqrt(qT)-F0)/sqrt(qT)
    XX2=(FF(x1,x2,Q2,tau2,qT,K,K2)/sqrt(qT+K2)-F0)/sqrt(qT-K2)
    XX3=(FF(x1,x2,Q2,tau2,qT,K,K3)/sqrt(qT+K3)-F0)/sqrt(qT-K3)
    XX4=(FF(x1,x2,Q2,tau2,qT,K,K4)/sqrt(qT+K4)-F0)/sqrt(qT-K4)
    XX5=(FF(x1,x2,Q2,tau2,qT,K,K5)/sqrt(qT+K5)-F0)/sqrt(qT-K5)
    XX6=(FF(x1,x2,Q2,tau2,qT,K,K6)/sqrt(qT+K6)-F0)/sqrt(qT-K6)
    XX7=(FF(x1,x2,Q2,tau2,qT,K,K7)/sqrt(qT+K7)-F0)/sqrt(qT-K7)
    XX8=(FF(x1,x2,Q2,tau2,qT,K,K8)/sqrt(qT+K8)-F0)/sqrt(qT-K8)
    XX9=0d0
    
    valueMax=deltaK*(XX1+4d0*XX2+2d0*XX3+4d0*XX4&
            +2d0*XX5+4d0*XX6+2d0*XX7+4d0*XX8+XX9)/24d0
    
    integraldK=valueMax+sqrt(4d0*qT)*F0
    
!     write(*,*) '>>>',XX1,XX2,XX3,XX4,XX5,'--->',valueMax,'-->',integraldK

    
    !!! k<0 integral
    deltaK=-qT
    K2=deltaK/4d0
    K3=deltaK/2d0
    K4=-qT-deltaK/4d0
    
    F0=FF(x1,x2,Q2,tau2,qT,K,-qT)/sqrt(2d0*qT)
    
    XX1=(FF(x1,x2,Q2,tau2,qT,K,0d0)/sqrt(qT)-F0)/sqrt(qT)
    XX2=(FF(x1,x2,Q2,tau2,qT,K,K2)/sqrt(qT-K2)-F0)/sqrt(qT+K2)
    XX3=(FF(x1,x2,Q2,tau2,qT,K,K3)/sqrt(qT-K3)-F0)/sqrt(qT+K3)
    XX4=(FF(x1,x2,Q2,tau2,qT,K,K4)/sqrt(qT-K4)-F0)/sqrt(qT+K4)
    XX5=(FF(x1,x2,Q2,tau2,qT,K,K5)/sqrt(qT-K5)-F0)/sqrt(qT+K5)
    XX6=(FF(x1,x2,Q2,tau2,qT,K,K6)/sqrt(qT-K6)-F0)/sqrt(qT+K6)
    XX7=(FF(x1,x2,Q2,tau2,qT,K,K7)/sqrt(qT-K7)-F0)/sqrt(qT+K7)
    XX8=(FF(x1,x2,Q2,tau2,qT,K,K8)/sqrt(qT-K8)-F0)/sqrt(qT+K8)
    XX9=0d0
        
!     write(*,*) '>>>',XX1,XX2,XX3,XX4,XX5,'--->',valueMax,'-->',integraldK
    
    valueMax=deltaK*(XX1+4d0*XX2+2d0*XX3+4d0*XX4&
            +2d0*XX5+4d0*XX6+2d0*XX7+4d0*XX8+XX9)/24d0
    
    
    integraldK=integraldK+(valueMax+sqrt(4d0*qT)*F0)

end function integraldK



!!! to integrate over KK I split the integral {qT,tau} and in each sector subtract the singularity
!!! the formula is 
!!! \int_{qT}^{tau} f(k)/sqrt(k-qT)/sqrt(tau-k) dk = 
!!!     =int_{qT}^qX (f(k)/sqrt(tau-k)-f(qT)/sqrt(tau-qT))/sqrt(k-qT)+2sqrt(qX-qT)/sqrt(tau-qT)f(qT)
!!!         +int_{qX}^tau (f(k)/sqrt(k-qT)-f(tau)/sqrt(tau-qT))/sqrt(tau-k)+2sqrt(tau-qX)/sqrt(tau-qT)f(tau)
function integralKK(x1,x2,Q2,tau2,qT)
    real*8::integralKK,x1,x2,Q2,tau2,qT,tau,qX
    real*8::deltaK,K2,K3,K4,K5,K6,K7,K8
    real*8::XX1,XX2,XX3,XX4,XX5,XX6,XX7,XX8,XX9,F0,valueMax
    
    tau=sqrt(tau2)
    qX=(qT+tau)/2d0
    
!     do i=0,20
!     write(*,*) "{",tau*i/20d0+(1d0-i/20d0)*qT,",",integraldK(x1,x2,Q2,tau2,qT,tau*i/20d0+(1d0-i/20d0)*qT),"},"
!     end do
!     write(*,*) 'qX', qX,qT,tau
    !!! k<qX integral
    deltaK=qX-qT
    K2=qT+deltaK/8d0
    K3=qT+deltaK/4d0
    K4=qT+3d0*deltaK/8d0
    K5=qT+4d0*deltaK/8d0
    K6=qT+5d0*deltaK/8d0
    K7=qT+6d0*deltaK/8d0
    K8=qT+7d0*deltaK/8d0
    
!     write(*,*) '-------',K2,K3,K4
    
    F0=integraldK(x1,x2,Q2,tau2,qT,qT)/sqrt(tau-qT)
    
    XX1=0d0
    XX2=(integraldK(x1,x2,Q2,tau2,qT,K2)/sqrt(tau-K2)-F0)/sqrt(K2-qT)
    XX3=(integraldK(x1,x2,Q2,tau2,qT,K3)/sqrt(tau-K3)-F0)/sqrt(K3-qT)
    XX4=(integraldK(x1,x2,Q2,tau2,qT,K4)/sqrt(tau-K4)-F0)/sqrt(K4-qT)
    XX5=(integraldK(x1,x2,Q2,tau2,qT,K5)/sqrt(tau-K5)-F0)/sqrt(K5-qT)
    XX6=(integraldK(x1,x2,Q2,tau2,qT,K6)/sqrt(tau-K6)-F0)/sqrt(K6-qT)
    XX7=(integraldK(x1,x2,Q2,tau2,qT,K7)/sqrt(tau-K7)-F0)/sqrt(K7-qT)
    XX8=(integraldK(x1,x2,Q2,tau2,qT,K8)/sqrt(tau-K8)-F0)/sqrt(K8-qT)
    XX9=(integraldK(x1,x2,Q2,tau2,qT,qX)/sqrt(tau-qX)-F0)/sqrt(qX-qT)
    
    valueMax=deltaK*(XX1+4d0*XX2+2d0*XX3+4d0*XX4&
            +2d0*XX5+4d0*XX6+2d0*XX7+4d0*XX8+XX9)/24d0
    
    integralKK=valueMax+2d0*sqrt(qX-qT)*F0
    
!     write(*,*) 'XXX',XX1,XX2,XX3,XX4,XX5,integralKK
!     write(*,*) '>>>',valueMax,2d0*sqrt(qX-qT)*F0
    
        !!! k>qX integral
    deltaK=tau-qX
    K2=qX+deltaK/8d0
    K3=qX+deltaK/4d0
    K4=qX+3d0*deltaK/8d0
    K5=qX+4d0*deltaK/8d0
    K6=qX+5d0*deltaK/8d0
    K7=qX+6d0*deltaK/8d0
    K8=qX+7d0*deltaK/8d0
    
    F0=integraldK(x1,x2,Q2,tau2,qT,tau)/sqrt(tau-qT)
    
    XX1=(integraldK(x1,x2,Q2,tau2,qT,qX)/sqrt(qX-qT)-F0)/sqrt(tau-qX)
    XX2=(integraldK(x1,x2,Q2,tau2,qT,K2)/sqrt(K2-qT)-F0)/sqrt(tau-K2)
    XX3=(integraldK(x1,x2,Q2,tau2,qT,K3)/sqrt(K3-qT)-F0)/sqrt(tau-K3)
    XX4=(integraldK(x1,x2,Q2,tau2,qT,K4)/sqrt(K4-qT)-F0)/sqrt(tau-K4)
    XX5=(integraldK(x1,x2,Q2,tau2,qT,K5)/sqrt(K5-qT)-F0)/sqrt(tau-K5)
    XX6=(integraldK(x1,x2,Q2,tau2,qT,K6)/sqrt(K6-qT)-F0)/sqrt(tau-K6)
    XX7=(integraldK(x1,x2,Q2,tau2,qT,K7)/sqrt(K7-qT)-F0)/sqrt(tau-K7)
    XX8=(integraldK(x1,x2,Q2,tau2,qT,K8)/sqrt(K8-qT)-F0)/sqrt(tau-K8)
    XX9=0d0
    
    valueMax=deltaK*(XX1+4d0*XX2+2d0*XX3+4d0*XX4&
            +2d0*XX5+4d0*XX6+2d0*XX7+4d0*XX8+XX9)/24d0
    
!     write(*,*) 'XXX',XX1,XX2,XX3,XX4,XX5,valueMax+2d0*sqrt(tau-qX)*F0
!     write(*,*) '>>>',valueMax,2d0*sqrt(tau-qX)*F0
    
    integralKK=integralKK+valueMax+2d0*sqrt(tau-qX)*F0

end function integralKK


!!!! I use the angle representtion of the integral
!!!! ch=Cosh[omega]
function integraldK_K15(x1,x2,Q2,tau2,qT,ch)
    real*8::integraldK_K15,x1,x2,Q2,tau2,qT,ch
    real*8::fff,cc,inter
    integer::i
    
    !delta=piHalf
    !av=piHalf
    
    inter=0d0
    do i=1,15
        cc=cos(Xi_k15(i)*piHalf+piHalf)
        fff=qT**2/2d0*(ch**2-cc**2)*FF(x1,x2,Q2,tau2,qT,qt*ch,qt*cc)
        inter=inter+Wi_k15(i)*fff
    end do
        
    integraldK_K15=piHalf*inter
    

end function integraldK_K15

!!!!! Qupper is the optional argument for the cut of the upper limit
function integralKK_GK(x1,x2,Q2,tau2,qT,Qupper,HasEndPointSingularity)
    real*8::integralKK_GK,x1,x2,Q2,tau2,qT
    real*8,optional::Qupper
    logical,optional::HasEndPointSingularity
    logical::hasSingularity
    real*8::ch,Q,upperL
    real*8,parameter::cutQfactor=5d0!!!!! I cut upper limit by cutQfactor*Q (5 gives ~1-2% accuracy)
    real*8,parameter::tolerance=0.001d0
    real*8::delta,av,g7,k15,fff,eps,cutPart,cutLimit
    integer::i
    
    if(present(HasEndPointSingularity)) then
        hasSingularity=HasEndPointSingularity
    else
        hasSingularity=.false.
    end if
    
    if(present(Qupper)) then
        Q=Qupper        
        upperL=log(Q/qT+sqrt((Q/qT)**2-1d0))
    else
        Q=sqrt(tau2)
        upperL=log(cutQfactor*Q/qT+sqrt((cutQfactor*Q/qT)**2-1d0))
    end if
    
    !!!! if the integral has HasEndPointSingularity then 
    !!!! I cut the last 10% of the integrla and integrate it with Gauss-Chebishev routine
    if(hasSingularity) then
        cutLimit=0.9d0*upperL
    else
        cutLimit=upperL
    end if
   
    delta=cutLimit/2d0
    av=cutLimit/2d0
    
    g7=0d0
    k15=0d0
    do i=1,15
        ch=cosh(Xi_k15(i)*delta+av)
        fff=integraldK_K15(x1,x2,Q2,tau2,qT,ch)
        g7=g7+Wi_g7(i)*fff
        k15=k15+Wi_k15(i)*fff        
    end do
    
    eps=delta*abs(k15)*tolerance
    
    if(hasSingularity) then
        cutPart=integralKK_GC(x1,x2,Q2,tau2,qT,qT*cosh(cutLimit),qT*cosh(upperL),eps)
    else
        cutPart=0d0
    end if
    
    if(delta*abs(k15-g7)>eps) then
        integralKK_GK=integralKK_GK_rec(x1,x2,Q2,tau2,qT,0d0,av,eps)+integralKK_GK_rec(x1,x2,Q2,tau2,qT,av,cutLimit,eps)+cutPart
    else
        integralKK_GK=delta*k15+cutPart
    end if
    !write(*,*) '------------>',integralKK_GK-cutPart,cutPart

end function integralKK_GK

recursive function integralKK_GK_rec(x1,x2,Q2,tau2,qT,w1,w2,eps) result(res)
    real*8::res,x1,x2,Q2,tau2,qT,w1,w2,eps
    real*8::delta,av,g7,k15,ch,fff
    integer::i
        
    delta=(w2-w1)/2d0
    av=(w2+w1)/2d0
    
    g7=0d0
    k15=0d0
    do i=1,15
        ch=cosh(Xi_k15(i)*delta+av)
        fff=integraldK_K15(x1,x2,Q2,tau2,qT,ch)
        g7=g7+Wi_g7(i)*fff
        k15=k15+Wi_k15(i)*fff
    end do
    
    if(delta*abs(k15-g7)>eps) then
        res=integralKK_GK_rec(x1,x2,Q2,tau2,qT,w1,av,eps)+integralKK_GK_rec(x1,x2,Q2,tau2,qT,av,w2,eps)
    else
        res=delta*k15
    end if

end function integralKK_GK_rec


!!!!! Gauss-Chebichev integral from K1, K2
recursive function integralKK_GC(x1,x2,Q2,tau2,qT,K1,K2,eps) result(res)
    real*8::res,x1,x2,Q2,tau2,qT,K1,K2
    real*8::delta,av,g7,k15,fff,eps,ch,cc,tau,Jacobian
    integer::i
        
    delta=piHalf
    av=piHalf
    tau=sqrt(tau2)
    
    g7=0d0
    k15=0d0
    do i=1,15
        cc=cos(Xi_k15(i)*delta+av)
        ch=(K1*(cc+1d0)+K2*(1d0-cc))/(2d0*qT)
        Jacobian=(K2-K1)*sqrt(1d0-cc**2)/sqrt((K1+K2+(K1-K2)*cc)**2-4d0*qT**2)        
        fff=Jacobian*integraldK_K15(x1,x2,Q2,tau2,qT,ch)
        g7=g7+Wi_g7(i)*fff
        k15=k15+Wi_k15(i)*fff        
    end do
        
    if(delta*abs(k15-g7)>eps) then
        res=integralKK_GC(x1,x2,Q2,tau2,qT,K1,(K1+K2)/2d0,eps)+integralKK_GC(x1,x2,Q2,tau2,qT,(K1+K2)/2d0,K2,eps)
    else
        res=delta*k15
    end if

end function integralKK_GC

end module ONE

!----------------------------------------------------------------------------------------

program example
use aTMDe_control
use aTMDe_Numerics
use ONE
implicit none

real*8,parameter::hc2=0.38937933800000002d0
real*8::Q2,qT2,tau,tau2,qT,s,y,x1,x2,temp
integer::i
integer,parameter::orderH_global=2

call artemide_Initialize('const-TMD-inKT_NNLO',prefix='Prog/PowerCorr/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.253434d0, 9.04351d0, 346.999d0, 2.47992d0, -5.69988d0, 0.d0, 0.d0/))

! Q2=91d0**2
! qT2=1d0**2
! s=1800d0**2
s=27.43d0**2
Q2=7d0**2
qT2=1d0
y=0d0
tau=sqrt(Q2+qT2)
tau2=Q2+qT2
x1=sqrt((Q2+qT2)/s)*exp(y)
x2=sqrt((Q2+qT2)/s)*exp(-y)

! call SetM1(2)
! 
! do i=0,40
!     qT=i/40d0*Sqrt(Q2)+0.0001d0
! 
!     qT2=qT**2
!     tau=sqrt(Q2+qT2)
!     tau2=Q2+qT2
!     x1=sqrt((Q2+qT2)/s)*exp(y)
!     x2=sqrt((Q2+qT2)/s)*exp(-y)
! 
!     temp=integralKK_GK(x1,x2,Q2,tau2,qT,Qupper=tau,HasEndPointSingularity=.true.)
!     write(*,'("{",F8.3,",",F14.12,"},")') qT,pi*temp
! end do


! write(*,*) temp

!     do i=-20,20
!     write(*,*) "{",qT*i/20d0,",",FF(x1,x2,Q2,tau**2,qT,5d0,qT*i/20d0),"},"
!     end do

! temp=FF(x1,x2,Q2,tau**2,qT,5d0,0.5d0)


write(*,*) '---------------- Original xSec -----------------------'
call SetM1(2)

do i=0,80
qT=i/80d0*Sqrt(Q2)/2d0+0.0001d0

qT2=qT**2
tau=sqrt(Q2+qT2)
tau2=Q2+qT2
x1=sqrt((Q2+qT2)/s)*exp(y)
x2=sqrt((Q2+qT2)/s)*exp(-y)

temp=integralKK_GK(x1,x2,Q2,tau2,qT,HasEndPointSingularity=.false.)
write(*,'("{",F8.3,",",F14.12,"},")') qT,pi*temp
end do

write(*,*) '---------------- Original xSec cut by tau-------------'
call SetM1(2)

do i=0,80
qT=i/80d0*Sqrt(Q2)/2d0+0.0001d0

qT2=qT**2
tau=sqrt(Q2+qT2)
tau2=Q2+qT2
x1=sqrt((Q2+qT2)/s)*exp(y)
x2=sqrt((Q2+qT2)/s)*exp(-y)

temp=integralKK_GK(x1,x2,Q2,tau2,qT,Qupper=tau,HasEndPointSingularity=.false.)
write(*,'("{",F8.3,",",F14.12,"},")') qT,pi*temp
end do

write(*,*) '---------------- Power-suppressed term -----------------------'

call SetM1(1)

do i=0,80

qT=i/80d0*Sqrt(Q2)/2d0+0.0001d0

qT2=qT**2
tau2=Q2+qT2
tau=sqrt(Q2+qT2)
x1=sqrt((Q2+qT2)/s)*exp(y)
x2=sqrt((Q2+qT2)/s)*exp(-y)

temp=integralKK_GK(x1,x2,Q2,tau2,qT,Qupper=tau,HasEndPointSingularity=.true.)
write(*,'("{",F8.3,",",F14.12,"},")') qT,pi*temp
end do

write(*,*) '---------------- Full calculation -----------------------'

call SetM1(3)

do i=0,80

qT=i/80d0*Sqrt(Q2)/2d0+0.0001d0

qT2=qT**2
tau=sqrt(Q2+qT2)
tau2=Q2+qT2
x1=sqrt((Q2+qT2)/s)*exp(y)
x2=sqrt((Q2+qT2)/s)*exp(-y)

temp=integralKK_GK(x1,x2,Q2,tau2,qT,Qupper=tau,HasEndPointSingularity=.true.)
write(*,'("{",F8.3,",",F14.12,"},")') qT,pi*temp
end do

end program

