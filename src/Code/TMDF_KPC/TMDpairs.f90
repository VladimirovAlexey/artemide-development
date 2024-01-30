!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! List of functions to substitute into the cross-section
!!! These functions contain product of TMD distributions and all what is neccesary to attach
!!!
function TMD_pair(Q2,x1,x2,k1,k2,mu,process)
real(dp),intent(in)::Q2,x1,x2,k1,k2,mu
integer,dimension(1:3),intent(in)::process
real(dp)::TMD_pair
real(dp),dimension(-5:5)::FA,FB,FAB
real(dp)::param
integer::h1,h2

!increment counter
GlobalCounter=GlobalCounter+1
LocalCounter=LocalCounter+1

h1=process(1)
h2=process(2)


SELECT CASE(process(3))
  !!!test cases
  CASE(0,10000,20000,30000)
    !param=0.1d0
    !TMD_pair=(param/(k1+param**2)**1.5)*(param/(k2+param**2)**1.5)/(pix4)
    param=2.5d0
    TMD_pair=(3*param*(2*param**2-3*k1)/(k1+param**2)**3.5)*(3*param*(2*param**2-3*k2)/(k2+param**2)**3.5)/(pix4)
  CASE(9999,19999,29999,39999)
    TMD_pair=(Exp(-0.2d0*k1)+1/(k1+2.))*(Exp(-0.2d0*k2)+1/(k2+2.))
  CASE(9998,19998,29998,39998)
    TMD_pair=(Exp(-0.2d0*k1)+1/(k1+2.))*(Exp(-0.2d0*k2)+1/(k2+2.))

  CASE (1) !pp->gamma
    ! e_q^2 *F_q(A)*F_qbar(B)
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,h1)
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,-h2) !!! -h2, to myltiply quarks by anti-quarks in FAB
     FAB=FA*FB

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

!--------------------------------------------------------------------------------
  CASE (2,3,4,5) !Delta^{GG'}z_{+l}z_{+f}f1f1
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,h1)
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,-h2)!!! -h2, to myltiply quarks by anti-quarks in FAB
     FAB=FA*FB

    TMD_pair=XTMD_pairZpZp(FAB,Q2)



  CASE DEFAULT
    write(*,*) ErrorString('undefined process: ',moduleName),process
    write(*,*) color('Evaluation stop',c_red_bold)
    stop
 END SELECT


end function TMD_pair


!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!!! Combination Delta^{GG'} z_{+l}z_{+f} FF
function XTMD_pairZpZp(FAB,Q2)
     real(dp)::XTMD_pairZpZp
     real(dp),intent(in)::Q2
    !!cross-seciton parameters
     real(dp),dimension(-5:5),intent(in):: FAB

     !!!parameters of Z boson coupling
!      real(dp),parameter:: paramU=0.40329064872689d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=2/3
!      real(dp),parameter:: paramD=0.51983027428079d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1/3
!      real(dp),parameter:: paramS=0.51983027428079d0  !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1/3
!      real(dp),parameter:: paramC=0.40329064872689d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=2/3
!      real(dp),parameter:: paramB=0.51983027428079d0  !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1/3
!      real(dp),parameter:: paramL=0.35358707798999d0  !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1

     !!!parameters of Z-gamma boson coupling
!      real(dp),parameter:: paramMIXU=0.1515661518957d0 !! e(T3-2e sW^2)/2sWcW for eq=+2/3,, T3=+1/2
!      real(dp),parameter:: paramMIXD=0.1367184036034d0 !! e(T3-2e sW^2)/2sWcW for eq=-1/3,T3=-1/2
!      real(dp),parameter:: paramMIXS=0.1367184036034d0  !! e(T3-2e sW^2)/2sWcW for eq=-1/3, T3=-1/2
!      real(dp),parameter:: paramMIXC=0.1515661518957d0 !! e(T3-2e sW^2)/2sWcW for eq=+2/3,, T3=+1/2
!      real(dp),parameter:: paramMIXB=0.1367184036034d0  !! e(T3-2e sW^2)/2sWcW for eq=-1/3, T3=-1/2
!      real(dp),parameter:: paramMIXL=0.0445432448766d0  !! e(T3-2e sW^2)/2sWcW for eq=-1, T3=-1/2

     XTMD_pairZpZp=&
     (&!gamma-part
      4d0/9d0*FAB(2)&
      +1d0/9d0*FAB(1)&
      +1d0/9d0*FAB(3)&
      +4d0/9d0*FAB(4)&
      +1d0/9d0*FAB(5)&
      +4d0/9d0*FAB(-2)&
      +1d0/9d0*FAB(-1)&
      +1d0/9d0*FAB(-3)&
      +4d0/9d0*FAB(-4)&
      +1d0/9d0*FAB(-5))&
     +&!gamma-Z interference
     paramMIXL*(&
      paramMIXU*FAB(2)&
      +paramMIXD*FAB(1)&
      +paramMIXS*FAB(3)&
      +paramMIXC*FAB(4)&
      +paramMIXB*FAB(5)&
      +paramMIXU*FAB(-2)&
      +paramMIXD*FAB(-1)&
      +paramMIXS*FAB(-3)&
      +paramMIXC*FAB(-4)&
      +paramMIXB*FAB(-5))*&
      2d0*Q2*(Q2-MZ2)/((Q2-MZ2)**2+GammaZ2*MZ2)&
     +&!ZZ-contributions
     paramL*(&
      paramU*FAB(2)&
      +paramD*FAB(1)&
      +paramS*FAB(3)&
      +paramC*FAB(4)&
      +paramB*FAB(5)&
      +paramU*FAB(-2)&
      +paramD*FAB(-1)&
      +paramS*FAB(-3)&
      +paramC*FAB(-4)&
      +paramB*FAB(-5))*&
      Q2*Q2/((Q2-MZ2)**2+GammaZ2*MZ2)

end function XTMD_pairZpZp
