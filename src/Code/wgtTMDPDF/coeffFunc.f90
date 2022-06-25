!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of wgtTMDPDF module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for the mathing coefficient--------------------------------------
!!!--------------The order is accumulative pertrubative order of coefficient =1,2,3 (LO,NLO,NNLO)---------
!!!-------------------------------------------------------------------------------------------------------


 

  !!! the function which contains the functions of parameterizations
  function parametrizationString(z)
  real(dp)::z
  real(dp),dimension(1:parametrizationLength)::parametrizationString
      !lz=Log(z)
      !l1z=Log(1d0-z)
      parametrizationString=(/1d0/)
  
  end function parametrizationString
  
  !!! the function which contains 
    !!! int_z^1 parameterization at values of z -> 1
    !!! it is used to estimate integration error at z~1
  function parametrizationStringAt1(z)
  real(dp)::z
  real(dp),dimension(1:parametrizationLength)::parametrizationStringAt1
     parametrizationStringAt1=(/0d0/)
  
  end function parametrizationStringAt1

  
  !!!!Each coefficient is split to delta, sing x->1, regular
    
  !!!!!coefficient function q<-q delta-part
  function C_q_q_delta(alpha,Nf,Lmu)
  real(dp)::C_q_q_delta,Nf,alpha,Lmu
  
 C_q_q_delta=0d0
!  
!   if(order_global>=2) then      
!       C_q_q_delta=C_q_q_delta+alpha*(-4d0/3d0*zeta2-4d0*Lmu)
!   end if
  end function C_q_q_delta
  
  !!!!!coefficient function g<-g delta-part
  function C_g_g_delta(alpha,Nf,Lmu)
  real(dp)::C_g_g_delta,Nf,alpha,Lmu
  
  C_g_g_delta=0d0
 
!   if(order_global>=2) then
!   !!!checked (20.11.17) AV!!!
!       C_g_g_delta=C_g_g_delta+alpha*(-3d0*zeta2-(11d0-2d0/3d0*Nf)*Lmu)
!   end if  
  end function C_g_g_delta
  
  !!!!!coefficient function q<-q singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
  subroutine Set_CoeffSing1_q_q(alpha,Nf,Lmu)
  real(dp)::Nf,alpha,LLambda,Lmu,s1,s2,s3
    
  s1=0d0!!!coeff 1/(1-x)
  s2=0d0!!!coeff log(1-x)/(1-x)
  s3=0d0!!!coeff log(1-x)^2/(1-x)
  

!   if(order_global>=1) then    
!     s1=s1+alpha*(-16d0/3d0)*Lmu
!   end if
  
  CoeffSing1_q_q=(/s1,s2,s3/)
  
  end subroutine Set_CoeffSing1_q_q
  
  !!!!!coefficient function g<-g singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
  subroutine Set_CoeffSing1_g_g(alpha,Nf,Lmu)
  real(dp)::Nf,alpha,Lmu,s1,s2,s3
    
  s1=0d0!!!coeff 1/(1-x)
  s2=0d0!!!coeff log(1-x)/(1-x)
  s3=0d0!!!coeff log(1-x)^2/(1-x)
 
!   if(order_global>=1) then  
!   !!!checked (20.11.17) AV!!!
!     s1=s1+alpha*(-12d0)*Lmu
!   end if
  
  CoeffSing1_g_g=(/s1,s2,s3/)
  
  end subroutine Set_CoeffSing1_g_g
  
  !!!!!coefficient function q<-q regular-part  
  subroutine Set_Coeff_q_q(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is 1, it is WW-part of worm-gear function
  Coeff_q_q=(/1d0/) !1
!   if(order_global>=1) then
!           
!     Coeff_q_q=Coeff_q_q+alpha*(/&0d0/) !
!     
! !  write(*,*) 'regularPart=', regularPart/x
!   end if
  end subroutine Set_Coeff_q_q
 
  !!!!!coefficient function q<-g regular-part  
  subroutine Set_Coeff_q_g(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_g=(/0d0/)
!   if(order_global>=1) then
!     Coeff_q_g=Coeff_q_g+alpha*(/0d0/)!
!   end if
  end subroutine Set_Coeff_q_g
  
    !!!!!coefficient function g<-q regular-part  
  subroutine Set_Coeff_g_q(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_q=(/0d0/)
!   if(order_global>=1) then
!     Coeff_g_q=Coeff_g_q+alpha*(/0d0/)!
!     
!   end if
  end subroutine Set_Coeff_g_q
  
      !!!!!coefficient function g<-g regular-part  
  subroutine Set_Coeff_g_g(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_g=(/1d0/)
  
  if(outputLevel>0) write(*,*) &
                WarningString('gluon part of worm-gear T function is not known. Set alike quark.',moduleName)
  
  end subroutine Set_Coeff_g_g

   !!!!!coefficient function q<-qb regular-part  
  subroutine Set_Coeff_q_qb(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qb=(/0d0/)!
  
  end subroutine Set_Coeff_q_qb
  
     !!!!!coefficient function q<-qp regular-part  
  subroutine Set_Coeff_q_qp(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qp=(/0d0/)
  end subroutine Set_Coeff_q_qp
  
  !!! This function has been used during debuging
 subroutine CheckCoefficient(as,Nf,Lmu,z)
 real(dp)::Lmu,as,z,Nf,lz,l1z
 real(dp), dimension(1:23)::func
 real(dp), dimension(1:2)::func1
 
!   lz=Log(z)
!   l1z=Log(1d0-z)
!   func=(/l1z,l1z**2,l1z**3,&
! 	    1d0/z,lz/z,&  
! 	    lz,lz**2,lz**3,&
! 	    1d0,z,z**2,& 
! 	    z/(1d0-z)*lz, z*lz, (z**2)*lz,& 
! 	    z/(1d0-z)*lz**2,z*lz**2,& 
! 	    (lz/(1d0-z)+1d0)*l1z,lz*l1z,z*lz*l1z,& 
! 	    (1d0-z)/z*l1z, (1d0-z)*l1z, ((1d0-z)**2)*l1z,(1d0-z)*l1z**2/)
!      
!   func1=(/1d0/(1d0-z),Log(1d0-z)/(1d0-z)/)
!  
 !!Q->Q
!   call Set_CoeffSing1_q_q(as,Nf,Lmu) 
!   call Set_Coeff_q_q(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_q_q*func)+SUM(CoeffSing1_q_q*func1)
  
!   !!Q->G
!   call Set_Coeff_q_g(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_q_g*func)

!   !!Q->Q'
!   call Set_Coeff_q_qp(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_q_qp*func)

  !!Q->Qbar
!   call Set_Coeff_q_qb(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_q_qb*func)
 
!  !! G->Q
!   call Set_Coeff_g_q(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_g_q*func)

!	!!G->G
  call Set_CoeffSing1_g_g(as,Nf,Lmu) 
  call Set_Coeff_g_g(as,Nf,Lmu)  
  !write(*,*) SUM(Coeff_g_g*func)+SUM(CoeffSing1_g_g*func1)
 end subroutine CheckCoefficient
