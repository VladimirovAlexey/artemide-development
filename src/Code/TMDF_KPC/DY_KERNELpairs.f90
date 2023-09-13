!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! List of functions to substitute into the cross-section
!!! These functions contain kernels that multiply the TMD distributions in the convolution
!!!
!!! Q2, qT2 external variables
!!! x1,x2 -- Bjorken x's
!!! xi1,xi2 -- collinear momenta of TMDs
!!! k1,k2 -- transverse momenta SQUARED of TMDs
!!! cT -- Cos theta of Delta vector
!!! cA -- Cos alpha of Delta vector = tau2 Delta2
function DY_KERNEL_pair(Q2,qT2,x1,x2,xi1,xi2,k1,k2,cT,cA,process)
real(dp),intent(in)::Q2,qT2,x1,x2,xi1,xi2,k1,k2,cT,CA
integer,intent(in)::process
real(dp)::DY_KERNEL_pair
real(dp)::k1k2


SELECT CASE(process)

  CASE (1) !1

    DY_KERNEL_pair=1._dp

  CASE (2) ! EW DY ~P1 f1f1

    DY_KERNEL_pair=-xi1/xi2*k2/Q2

  CASE (3) ! EW DY ~P2 f1f1

    DY_KERNEL_pair=-xi2/xi1*k1/Q2

  CASE (4) ! EW DY ~P3 f1f1
    k1k2=(qT2-cA**2**(Q2-qT2*cT**2))/4
    DY_KERNEL_pair=2*k1k2/Q2

  CASE (5) ! EW DY ~P1 h1h1
    k1k2=(qT2-cA**2**(Q2-qT2*cT**2))/4
    DY_KERNEL_pair=xi1/xi2*k2/Q2*k1k2/M2

  CASE (6) ! EW DY ~P2 h1h1
    k1k2=(qT2-cA**2**(Q2-qT2*cT**2))/4
    DY_KERNEL_pair=xi2/xi1*k1/Q2*k1k2/M2

  CASE (7) ! EW DY ~P3 h1h1
    DY_KERNEL_pair=2*k1*k2/Q2/M2

  CASE (8) ! EW DY ~P1A f1f1
    DY_KERNEL_pair=xi1/x1-x2/xi2*k2/(Q2+qT2)

  CASE (9) ! EW DY ~P2A f1f1
    DY_KERNEL_pair=xi2/x2-x1/xi1*k1/(Q2+qT2)

  CASE DEFAULT
    write(*,*) ErrorString('undefined process 2 variables: ',moduleName),process
    write(*,*) color('Evaluation stop',c_red_bold)
    stop
 END SELECT


end function DY_KERNEL_pair
