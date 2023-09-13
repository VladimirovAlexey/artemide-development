!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! List of functions to substitute into the cross-section
!!! These functions contain product of TMD distributions and all what is neccesary to attach
!!!
function TMD_pair(Q2,x1,x2,k1,k2,mu,process)
real(dp),intent(in)::Q2,x1,x2,k1,k2,mu
integer,dimension(1:3),intent(in)::process
real(dp)::TMD_pair
real(dp),dimension(-5:5)::FA,FB,FAB

!increment counter
GlobalCounter=GlobalCounter+1
LocalCounter=LocalCounter+1

SELECT CASE(process(1))

  CASE (1) !pp->gamma
    ! e_q^2 *F_q(A)*F_qbar(B)
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(3))
     FAB=FA*(FB(5:-5:-1))

    TMD_pair=FAB(1)/9.d0&
      +FAB(2)*4.d0/9.d0&
      +FAB(3)/9.d0&
      +FAB(4)*4d0/9.d0&
      +FAB(5)/9d0&
      +FAB(-1)/9.d0&
      +FAB(-2)*4.d0/9.d0&
      +FAB(-3)/9.d0&
      +FAB(-4)*4d0/9.d0&
      +FAB(-5)/9d0
  CASE DEFAULT
    write(*,*) ErrorString('undefined process: ',moduleName),process
    write(*,*) color('Evaluation stop',c_red_bold)
    stop
 END SELECT


end function TMD_pair
