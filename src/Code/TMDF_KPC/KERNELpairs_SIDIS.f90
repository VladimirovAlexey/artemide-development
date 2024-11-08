!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! List of functions to substitute into the cross-section
!!! These functions contain kernels that multiply the TMD distributions in the convolution
!!!
!!! Q2, tau2, qT2 external variables
!!! S, Lambda -- integration variables
!!! sA=sin(alpha)
!!! cT=cos(theta)
!!! process -- integer that specifies process
function SIDIS_KERNEL(Q2,tau2,qT2,S,Lam,sA,cT,process)
real(dp),intent(in)::Q2,tau2,qT2,S,Lam,sA,cT
integer,intent(in)::process
real(dp)::SIDIS_KERNEL


SELECT CASE(process)

  CASE (2001) !unpol. DY, EW-procces

    SIDIS_KERNEL=1._dp

  CASE DEFAULT
    write(*,*) ErrorString('undefined process 2 variables: ',moduleName),process
    write(*,*) color('Evaluation stop',c_red_bold)
    stop
 END SELECT


end function SIDIS_KERNEL
