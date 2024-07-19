!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	This file contains the part of the code, which is common for all TMD-OPE modules
!       that operates at twist-2. It is inclucded (as a text).
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!	This part collect the definition of expressions used in the large-X resummation
!
!	v.3.00 Created (AV, 19.07.2024)
!
!				A.Vladimirov (19.07.2024)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! The large-X resummation leds to the following expression
!!! (\delta(1-x)-A/(1-x)^A) Exp[EXP]
!!! where A and EXP are perturbative expressions
!!! in this module create these expressions from the known part

!!!!!Argument of the exponent
!!! QUARK channel
pure function LargeX_EXP_q_q(alpha,Nf,Lmu)
  real(dp),intent(in)::alpha,Nf,Lmu
  real(dp)::LargeX_EXP_q_q
  real(dp)::dd1,dd2

  LargeX_EXP_q_q=0.d0

  !!! the coefficients of delta-function are extracted from coefficient files
  if(orderLX>=1) then
    dd1=C_q_q_delta_1(Nf,Lmu)
    LargeX_EXP_q_q=LargeX_EXP_q_q+alpha*dd1

  if(orderLX>=2) then
    dd2=C_q_q_delta_2(Nf,Lmu)
    LargeX_EXP_q_q=LargeX_EXP_q_q+alpha*alpha*(2*dd2-dd1**2)/2.d0

  if(orderLX>=3) then
    LargeX_EXP_q_q=LargeX_EXP_q_q+alpha**3*(3*C_q_q_delta_3(Nf,Lmu)-3*dd2*dd1+dd1**3)/3.d0
  end if
  end if
  end if
end function LargeX_EXP_q_q

!!! GLUON channel
pure function LargeX_EXP_g_g(alpha,Nf,Lmu)
  real(dp),intent(in)::alpha,Nf,Lmu
  real(dp)::LargeX_EXP_g_g
  real(dp)::dd1,dd2

  LargeX_EXP_g_g=0.d0

  !!! the coefficients of delta-function are extracted from coefficient files
  if(orderLX>=1) then
    dd1=C_g_g_delta_1(Nf,Lmu)
    LargeX_EXP_g_g=LargeX_EXP_g_g+alpha*dd1

  if(orderLX>=2) then
    dd2=C_g_g_delta_2(Nf,Lmu)
    LargeX_EXP_g_g=LargeX_EXP_g_g+alpha*alpha*(2*dd2-dd1**2)/2.d0

  if(orderLX>=3) then
    LargeX_EXP_g_g=LargeX_EXP_g_g+alpha**3*(3*C_g_g_delta_3(Nf,Lmu)-3*dd2*dd1+dd1**3)/3.d0
  end if
  end if
  end if
end function LargeX_EXP_g_g
