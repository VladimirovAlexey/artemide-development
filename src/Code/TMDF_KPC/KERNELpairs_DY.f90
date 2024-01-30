!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! List of functions to substitute into the cross-section
!!! These functions contain kernels that multiply the TMD distributions in the convolution
!!!
!!! Q2, tau2, qT2 external variables
!!! S, Lambda -- integration variables
!!! process -- integer that specifies process
function DY_KERNEL(Q2,tau2,qT2,S,Lam,process)
real(dp),intent(in)::Q2,tau2,qT2,S,Lam
integer,intent(in)::process
real(dp)::DY_KERNEL


SELECT CASE(process)

  CASE (1,2) !unpol. DY, EW-procces

    DY_KERNEL=1._dp

  CASE (3) !Z-boson ~P1 f1f1

    DY_KERNEL=tau2*(-1+S**2+Lam)/Q2/2

  CASE (4) !Z-boson ~P2 f1f1

    DY_KERNEL=tau2*Sqrt(Lam)*S/Q2

  CASE (5) !Z-boson ~P3 f1f1

    DY_KERNEL=1+tau2*(-1+S**2+Lam)/Q2/2

  CASE (6) !Z-boson ~P1A {f1f1}_A

    DY_KERNEL=-2*S*sqrt(tau2/Q2)

  CASE (7) !Z-boson ~P2A {f1f1}_A

    DY_KERNEL=-2*sqrt(Lam*tau2/Q2)

  CASE (8) !Z-boson ~P1 h1h1

    DY_KERNEL=tau2*(-1+S**2+Lam)*(2*Q2+tau2*(-1+S**2-Lam))/8/Q2/M2

  CASE (9) !Z-boson ~P2 h1h1

    DY_KERNEL=tau2*S*sqrt(Lam)*(2*Q2+tau2*(-1+S**2-Lam))/4/Q2/M2

  CASE (10) !Z-boson ~P3 h1h1

    DY_KERNEL=tau2**2*((S-1)**2-Lam)*((S+1)**2-Lam)/8/Q2/M2

  CASE (11) !Z-boson ~P4 {h1h1}_A

    DY_KERNEL=-tau2**(1.5)*Sqrt(Lam)*(Q2*(1+S**2-Lam)+tau2*(-1+S**2+Lam))/4/Q2/M2/sqrt(qT2)

  CASE (12) !Z-boson ~P5 {h1h1}_A

    DY_KERNEL=-tau2**(1.5)*S*(Q2*(-1+S**2-Lam)+2*Lam*tau2)/4/Q2/M2/sqrt(qT2)

  CASE DEFAULT
    write(*,*) ErrorString('undefined process 2 variables: ',moduleName),process
    write(*,*) color('Evaluation stop',c_red_bold)
    stop
 END SELECT


end function DY_KERNEL
