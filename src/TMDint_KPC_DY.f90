!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.00
!
!    Evaluation of the convolution integral for DY TMD cross-section with KPC
!
!    if you use this module please, quote 2307.13054
!
!    ver 3.0: created (AV, 07.09.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDint_KPC_DY
use aTMDe_Numerics
use IntegrationRoutines
use IO_functions
use TMDs_inKT

implicit none

private

character (len=12),parameter :: moduleName="TMDin-KPC-DY"
character (len=5),parameter :: version="v3.00"
!Last appropriate verion of constants-file
integer,parameter::inputver=30

public::KPC_DYconv

!increment counters
integer::GlobalCounter=0 !!!total counter of calls of TMD pairs
integer::LocalCounter=0 !!!counter of calls of TMD pairs within the current integrand

contains

INCLUDE 'Code/TMDint_KPC/DYTMDpairs.f90'

!!!--------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Function that computes the integral for KPC convolution in DY
!!! proc1 = (int,int,int) is the process def for TMD*TMD
!!! proc2 = int is the process definition for the integral kernel
!!! Q, qT, x1,x2, mu are usual DY variables
!!! THIS IS A SYMMETRIC VERSION (i.e. it should contain only cos(theta)
function KPC_DYconv(Q,qT,x1,x2,mu,proc1,proc2)
real(dp),intent(in)::Q,qT,x1,x2,mu
integer,intent(in),dimension(1:3)::proc1
integer,intent(in)::proc2
real(dp)::KPC_DYconv

real(dp)::tau2,deltaT

LocalCounter=0

tau2=Q**2+qT**2
deltaT=qT**2/tau2

KPC_DYconv=Integrate_GK(InTheta,0._dp,pi,0.0001_dp)
write(*,*) "LC=",LocalCounter

contains

function InTheta(theta)
real(dp)::InTheta
real(dp),intent(in)::theta
real(dp)::cT

real(dp)::S1,xi11,xi21,K11,K21

cT=cos(theta)

!!!! these are values of the parameters at deltaD=1
!!!! they can be passed further in order not to repeat the computation
S1=Sqrt(deltaT)*cT
xi11=x1/2*(1+S1)
xi21=x2/2*(1-S1)
K11=tau2/4*(1+S1)**2
K21=tau2/4*(1-S1)**2

InTheta=INT_overDelta(Q,tau2,deltaT,x1,x2,mu,proc1,proc2,cT,xi11,xi21,K11,K21)&
            +DY_TMD_pair(Q,xi11,xi21,k11,k21,mu,proc1) !! subtraction term
end function InTheta

end function KPC_DYconv

function INT_overDelta(Q,tau2,deltaT,x1,x2,mu,proc1,proc2,cT,xi11,xi21,K11,K21)
real(dp),intent(in)::Q,tau2,deltaT,x1,x2,mu,cT,xi11,xi21,K11,K21
integer,intent(in),dimension(1:3)::proc1
integer,intent(in)::proc2
real(dp)::INT_overDelta

INT_overDelta=Integrate_GK(InDelta,0._dp,1._dp,0.0001_dp)

contains

function InDelta(deltaD)
real(dp)::InDelta
real(dp),intent(in)::deltaD
real(dp)::S,Lam,xi1,xi2,K1,K2

S=Sqrt(deltaT*deltaD)*cT
Lam=(1-deltaT)*(1-deltaD)

xi1=x1/2*(1+S+sqrt(Lam))
xi2=x2/2*(1-S+sqrt(Lam))
K1=tau2/4*((1+S)**2-Lam)
K2=tau2/4*((1-S)**2-Lam)

!!! note that the values of parameters at d=1, computed in the earlier integral
!InDelta=(F1(xi1,K1)*F2(xi2,K2)-F1(xi11,K11)*F2(xi21,K21))/sqrt(1-deltaD)/2
InDelta=(DY_TMD_pair(Q,xi1,xi2,k1,k2,mu,proc1)-DY_TMD_pair(Q,xi11,xi21,k11,k21,mu,proc1))/sqrt(1-deltaD)/2


end function InDelta

end function INT_overDelta

!!!--------------------------------------------------------------------------------------------------

end module TMDint_KPC_DY
