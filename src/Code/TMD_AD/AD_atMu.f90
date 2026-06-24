!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of TMD_AD module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Anomalous dimensions at scale mu------------------------------------------
!!!-------------------------------------------------------------------------------------------------------
  
! The TMD rapidity anomalous dimension D, perturbative expression
! order=0 Dpert=as^0=0
! order=1 Dpert=as^1=0
! ...
!!!!! if order is not specified the orderD is taken
function Dpert(mu,bT,f)
    real(dp),intent(in):: bT,mu
    integer,intent(in):: f
    real(dp)::Dpert
    integer:: Nf
    real(dp)::LL,astrong

    LL=2_dp*LOG(bt*mu*C0_inv_const)
    astrong=As(mu)
    !Dpert=0_dp
    
    Nf=ActiveNf(mu)

    Dpert=Dpert_atL(astrong,Nf,LL,orderD,f)
end function Dpert

!!!!!! The TMD rapidity anomalous dimension D, perturbative expression
!!!!!! computed at (alpha,Nf,Lmu) and the order orderIN
pure function Dpert_atL(alpha,Nf,Lmu,orderIN,f)
    real(dp),intent(in):: alpha,Lmu
    integer,intent(in):: f,Nf,orderIN
    real(dp)::Dpert_atL
    integer:: n,k
    real(dp)::inter

    !!!! Dpert_atL = sum_{n=1}^order sum_{k=0}^n as^n L^k d_{nk}
    !!!! note that dnk =0 for k>n, so the double-sum is triangular
    !!!! the summation is made by Horner method such that it is more efficient.
    !!!! Horner method is a i1+a^2 i2+a^3 i3+.... = a(i1+a(i2+a(i3+...))), it saves multiplications dramatically
    Dpert_atL=0._dp
    if(f==0) then!gluon
      do n=orderIN,1,-1
        inter=0._dp
        do k=n,1,-1
            inter=(dnk_G(n,k,Nf)+inter)*Lmu
        end do
        inter=inter+dnk_G(n,0,Nf)

        Dpert_atL=alpha*(Dpert_atL+inter)
      end do
     else!quark
      do n=orderIN,1,-1
        inter=0._dp
        do k=n,1,-1
            inter=(dnk_Q(n,k,Nf)+inter)*Lmu
        end do
        inter=inter+dnk_Q(n,0,Nf)

        Dpert_atL=alpha*(Dpert_atL+inter)
      end do
    end if

end function Dpert_atL
  
! The TMD UV anomalous dimension gammaV
! order=0 gammaV=as^0=0
! order=1 gammaV=as^1
! ....
function gammaV(mu,f)
    real(dp),intent(in)::mu
    integer,intent(in)::f
    real(dp):: gammaV
    integer:: Nf,i    
    real(dp) ::astrong
    
    astrong=As(mu)
    
    Nf=ActiveNf(mu)
   
    gammaV=0_dp
    if(f==0) then   !!gluon case
     do i=orderV,1,-1
        gammaV=(gammaV+gammaV_G(i,Nf))*astrong
     end do
    else  !! quark case
     do i=orderV,1,-1
        gammaV=(gammaV+gammaV_Q(i,Nf))*astrong
     end do
    end if

end function gammaV

! The Gamma Cusp anomalous dimension
! order=0 gammaCusp=as^1
! order=1 gammaCusp=as^2
! ....
function gammaCUSP(mu,f)
    real(dp),intent(in)::mu
    integer,intent(in)::f
    real(dp) :: gammaCUSP
    integer:: Nf,i
    real(dp) ::astrong
    
    astrong=As(mu)
    
    Nf=activeNf(mu)
    
    gammaCUSP=0_dp
    if(f==0) then   !!gluon case
     do i=orderCusp,0,-1
        gammaCUSP=(gammaCUSP+GammaCusp_G(i,Nf))*astrong
     end do
    else  !! quark case
     do i=orderCusp,0,-1
        gammaCUSP=(gammaCUSP+GammaCusp_Q(i,Nf))*astrong
     end do
    end if
   
end function gammaCUSP

!!the resummed version of rapidity anomalous dimension
function Dresum(mu,bT,f)
    real(dp),intent(in)::mu,bT
    integer,intent(in)::f
    real(dp)::Dresum
    integer::Nf
    integer::n,k,l
    real(dp):: X,alpha,lX,a1X,LL,commulant,inter
    real(dp)::Xpow(0:orderDresum),lXpow(0:orderDresum)
    
    alpha=As(mu)
    Nf=ActiveNf(mu)
    LL=2_dp*LOG(bT*mu*C0_inv_const)
    
    X=betaQCD(0,Nf)*alpha*LL
    
    
    !!! in many models D evaluated at X=0, to speed up these computations I explicitely state it
    if(abs(X)<1d-7) then
    !!!!! This is optimized sum
    !!!! Dresum = a^2(d2+a*(d3+a*(d4+...)))=a^2 d2+a^3 d3+a^4 d4+....

    Dresum=0_dp
    if(f==0) then !! gluon case
        do n=orderDresum,2,-1
            Dresum=alpha*(dnk_G(n,0,Nf)+Dresum)
        end do
        Dresum=alpha*Dresum
    else !! quark case
        do n=orderDresum,2,-1
            Dresum=alpha*(dnk_Q(n,0,Nf)+Dresum)
        end do
        Dresum=alpha*Dresum
    end if
    
    !!! complete case
    else    
    
    !!!! pre compute powers
    lX=Log(1_dp-X)
    a1X=alpha/(1_dp-X)
    Xpow(0)=1d0
    lXpow(0)=1d0
    do k=1,orderDresum
        Xpow(k)=Xpow(k-1)*X
        lXpow(k)=lXpow(k-1)*lX
    end do
    
    !!!!! This is optimized sum
    !!!! commulant = lX+ a1X inter1+a1X^2 inter2+a1X^3 inter3+...
    if(f==0) then !! gluon case
        commulant=0._dp
        do n=orderDresum,1,-1
        inter=0_dp
            do k=0,n
            do l=0,n
                inter=inter+Xpow(k)*lXpow(l)*dnkl_G(n,k,l,Nf)
            end do
            end do
        commulant=a1X*(commulant+inter)
        end do
        commulant=commulant+lX
        
        !!!! thanks to Patricia Gutierrez for spoting error in this formula!!
        Dresum=-GammaCusp_G(0,Nf)/betaQCD(0,Nf)/2_dp*commulant
    else !!! quark case
        commulant=0._dp
        do n=orderDresum,1,-1
        inter=0_dp
            do k=0,n
            do l=0,n
                inter=inter+Xpow(k)*lXpow(l)*dnkl_Q(n,k,l,Nf)
            end do
            end do
        commulant=a1X*(commulant+inter)
        end do
        commulant=commulant+lX
        
        Dresum=-GammaCusp_Q(0,Nf)/betaQCD(0,Nf)/2_dp*commulant  
    end if
    
    if(ISNAN(Dresum)) then
        write(*,*) ErrorString('Dresum is NaN.',moduleName)
        write(*,*) 'At mu=',mu,'b=',bT,'Lmu=',2_dp*Log(mu*bT*C0_inv_const), 'X=',X,'log(1-x)=',lX
        write(*,*) 'Evaluation STOP'
        stop
    end if
    
    end if
end function Dresum

!-------------------zeta-lines -------------------------------------
!the exact zeta-line is in TMDR (because it depends on DNP)

!! the value of zeta_mu in the pertrubation theory with ri=0
function zetaMUpert(mu,bt,f)
  real(dp),intent(in)::mu,bT
  integer,intent(in)::f
  real(dp):: zetaMUpert
  integer::Nf,n,k
  real(dp)::alpha,LL,val,iter
  
  if(orderZETA>0) then !!! NLO,NNLO,...
  
    LL=2_dp*LOG(bt*mu*C0_inv_const)
    alpha=As(mu)
    Nf=ActiveNf(mu)    
    
    !!!! val = sum_{n=0}^order sum_{k=0}^(n+1) as^n L^k v_{nk}
    !!!! the summation is made by Horner method such that it is more efficient.
    !!!! Horner method is a i1+a^2 i2+a^3 i3+.... = a(i1+a(i2+a(i3+...))), it saves multiplications dramatically

    !!!!! Important!!!!
    !! the perturbative value for zeta, must be taken 1-order higher
    !! because there are double logarithms ~~beta0 L
    !! For 4-loop it requires 5-loop rad... The value which contains it is set to zero v(4,0)=0.
    if(f==0) then !!! gluon
        val=0._dp
        do n=orderZETA,1,-1
            iter = 0_dp
            do k=n+1,1,-1
                iter=(vnk_g(n,k,Nf)+iter)*LL
            end do
            iter=iter+vnk_g(n,0,Nf)
            val=alpha*(val+iter)
        end do   
        val=val+vnk_g(0,0,Nf)

    else !!!! quark
        val=0._dp
        do n=orderZETA,1,-1
            iter = 0_dp
            do k=n+1,1,-1
                iter=(vnk_q(n,k,Nf)+iter)*LL
            end do
            iter=iter+vnk_q(n,0,Nf)
            val=alpha*(val+iter)
        end do
        val=val+vnk_q(0,0,Nf)
    end if
  
    zetaMUpert=mu*C0_const/bT*EXP(-val)
  
  else if (orderZETA==0) then   !!! LO 
    zetaMUpert=mu*C0_const/bT
    
  else !!!! just in case
    zetaMUpert=1_dp
  end if
  
end function zetaMUpert
 
