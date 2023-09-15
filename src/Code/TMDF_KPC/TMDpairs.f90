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

!increment counter
GlobalCounter=GlobalCounter+1
LocalCounter=LocalCounter+1

SELECT CASE(process(1))
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

  CASE (2) !ppbar->gamma
    ! e_q^2 *F_q(A)*F_q(B)
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))
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
  CASE (3) !pp->Z
      !((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) *F_q(A)*F_qbar(B)
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))
     FAB=FA*(FB(5:-5:-1))

    TMD_pair=&
      FAB(1)*paramD&
      +FAB(2)*paramU&
      +FAB(3)*paramS&
      +FAB(4)*paramC&
      +FAB(5)*paramB&
      +FAB(-1)*paramD&
      +FAB(-2)*paramU&
      +FAB(-3)*paramS&
      +FAB(-4)*paramC&
      +FAB(-5)*paramB

!--------------------------------------------------------------------------------
  CASE (4) !ppbar->Z
      !((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) *F_q(A)*F_qbar(B)
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))
     FAB=FA*FB

    TMD_pair=&
      FAB(1)*paramD&
      +FAB(2)*paramU&
      +FAB(3)*paramS&
      +FAB(4)*paramC&
      +FAB(5)*paramB&
      +FAB(-1)*paramD&
      +FAB(-2)*paramU&
      +FAB(-3)*paramS&
      +FAB(-4)*paramC&
      +FAB(-5)*paramB

!--------------------------------------------------------------------------------
  CASE (5) !pp->Z+gamma
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))
     FAB=FA*(FB(5:-5:-1))

    TMD_pair=XTMD_pairForDYwithZgamma(FAB,Q2)
!--------------------------------------------------------------------------------
  CASE (6) !ppbar->Z+gamma
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))
     FAB=FA*FB

    !! we invert the order of FB
    TMD_pair=XTMD_pairForDYwithZgamma(FAB,Q2)
!--------------------------------------------------------------------------------
  CASE (7) !pp-> W+
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))

    TMD_pair=paramW_L*(&
    paramW_UD*(FA(2)*FB(-1)+FA(-1)*FB(2))&        !u*dbar+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(-3)*FB(2))&        !u*sbar+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(-5)*FB(2))&        !u*bbar+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(-1)*FB(4))&        !c*dbar+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(-3)*FB(4))&        !c*sbar+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(-5)*FB(4))&        !c*bbar+bbar*c
    )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------
  CASE (8) !pp-> W-
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))

    TMD_pair=paramW_L*(&
    paramW_UD*(FA(1)*FB(-2)+FA(-2)*FB(1))&        !d*ubar+ubar*d
    +paramW_US*(FA(3)*FB(-2)+FA(-2)*FB(3))&        !s*ubar+ubar*s
    +paramW_UB*(FA(5)*FB(-2)+FA(-2)*FB(5))&        !b*ubar+ubar*b
    +paramW_CD*(FA(1)*FB(-4)+FA(-4)*FB(1))&        !d*cbar+cbar*d
    +paramW_CS*(FA(3)*FB(-4)+FA(-4)*FB(3))&        !s*cbar+cbar*s
    +paramW_CB*(FA(5)*FB(-4)+FA(-4)*FB(5))&        !b*cbar+cbar*b
    )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------
  CASE (9) !pp-> W+ + W-
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))

    TMD_pair=paramW_L*(&
    paramW_UD*(FA(2)*FB(-1)+FA(1)*FB(-2)+FA(-2)*FB(1)+FA(-1)*FB(2))&    !u*dbar+d*ubar+ubar*d+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(3)*FB(-2)+FA(-2)*FB(3)+FA(-3)*FB(2))&    !u*sbar+s*ubar+ubar*s+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(5)*FB(-2)+FA(-2)*FB(5)+FA(-5)*FB(2))&    !u*bbar+b*ubar+ubar*b+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(1)*FB(-4)+FA(-4)*FB(1)+FA(-1)*FB(4))&    !c*dbar+d*cbar+cbar*d+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(3)*FB(-4)+FA(-4)*FB(3)+FA(-3)*FB(4))&    !c*sbar+s*cbar+cbar*s+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(5)*FB(-4)+FA(-4)*FB(5)+FA(-5)*FB(4))&    !c*bbar+b*cbar+cbar*b+bbar*c
    )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------
  CASE (10) !ppbar-> W+
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))
     FB=FB(5:-5:-1) !! inverse the quark order

    TMD_pair=paramW_L*(&
    paramW_UD*(FA(2)*FB(-1)+FA(-1)*FB(2))&        !u*dbar+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(-3)*FB(2))&        !u*sbar+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(-5)*FB(2))&        !u*bbar+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(-1)*FB(4))&        !c*dbar+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(-3)*FB(4))&        !c*sbar+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(-5)*FB(4))&        !c*bbar+bbar*c
    )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------
  CASE (11) !ppbar-> W-
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))
     FB=FB(5:-5:-1) !! inverse the quark order

    TMD_pair=paramW_L*(&
    paramW_UD*(FA(1)*FB(-2)+FA(-2)*FB(1))&        !d*ubar+ubar*d
    +paramW_US*(FA(3)*FB(-2)+FA(-2)*FB(3))&        !s*ubar+ubar*s
    +paramW_UB*(FA(5)*FB(-2)+FA(-2)*FB(5))&        !b*ubar+ubar*b
    +paramW_CD*(FA(1)*FB(-4)+FA(-4)*FB(1))&        !d*cbar+cbar*d
    +paramW_CS*(FA(3)*FB(-4)+FA(-4)*FB(3))&        !s*cbar+cbar*s
    +paramW_CB*(FA(5)*FB(-4)+FA(-4)*FB(5))&        !b*cbar+cbar*b
    )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------
  CASE (12) !ppbar-> W+ + W-
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))
     FB=FB(5:-5:-1) !! inverse the quark order

    TMD_pair=paramW_L*(&
    paramW_UD*(FA(2)*FB(-1)+FA(1)*FB(-2)+FA(-2)*FB(1)+FA(-1)*FB(2))&    !u*dbar+d*ubar+ubar*d+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(3)*FB(-2)+FA(-2)*FB(3)+FA(-3)*FB(2))&    !u*sbar+s*ubar+ubar*s+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(5)*FB(-2)+FA(-2)*FB(5)+FA(-5)*FB(2))&    !u*bbar+b*ubar+ubar*b+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(1)*FB(-4)+FA(-4)*FB(1)+FA(-1)*FB(4))&    !c*dbar+d*cbar+cbar*d+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(3)*FB(-4)+FA(-4)*FB(3)+FA(-3)*FB(4))&    !c*sbar+s*cbar+cbar*s+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(5)*FB(-4)+FA(-4)*FB(5)+FA(-5)*FB(4))&    !c*bbar+b*cbar+cbar*b+bbar*c
    )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)

!--------------------------------------------------------------------------------
  CASE (13) !pp-> W+
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))

    TMD_pair=&
    paramW_UD*(FA(2)*FB(-1)+FA(-1)*FB(2))&        !u*dbar+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(-3)*FB(2))&        !u*sbar+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(-5)*FB(2))&        !u*bbar+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(-1)*FB(4))&        !c*dbar+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(-3)*FB(4))&        !c*sbar+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(-5)*FB(4))        !c*bbar+bbar*c

!--------------------------------------------------------------------------------
  CASE (14) !pp-> W-
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))

    TMD_pair=&
    paramW_UD*(FA(1)*FB(-2)+FA(-2)*FB(1))&        !d*ubar+ubar*d
    +paramW_US*(FA(3)*FB(-2)+FA(-2)*FB(3))&        !s*ubar+ubar*s
    +paramW_UB*(FA(5)*FB(-2)+FA(-2)*FB(5))&        !b*ubar+ubar*b
    +paramW_CD*(FA(1)*FB(-4)+FA(-4)*FB(1))&        !d*cbar+cbar*d
    +paramW_CS*(FA(3)*FB(-4)+FA(-4)*FB(3))&        !s*cbar+cbar*s
    +paramW_CB*(FA(5)*FB(-4)+FA(-4)*FB(5))        !b*cbar+cbar*b

!--------------------------------------------------------------------------------
  CASE (15) !pp-> W+ + W-
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))

    TMD_pair=&
    paramW_UD*(FA(2)*FB(-1)+FA(1)*FB(-2)+FA(-2)*FB(1)+FA(-1)*FB(2))&    !u*dbar+d*ubar+ubar*d+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(3)*FB(-2)+FA(-2)*FB(3)+FA(-3)*FB(2))&    !u*sbar+s*ubar+ubar*s+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(5)*FB(-2)+FA(-2)*FB(5)+FA(-5)*FB(2))&    !u*bbar+b*ubar+ubar*b+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(1)*FB(-4)+FA(-4)*FB(1)+FA(-1)*FB(4))&    !c*dbar+d*cbar+cbar*d+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(3)*FB(-4)+FA(-4)*FB(3)+FA(-3)*FB(4))&    !c*sbar+s*cbar+cbar*s+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(5)*FB(-4)+FA(-4)*FB(5)+FA(-5)*FB(4))    !c*bbar+b*cbar+cbar*b+bbar*c

!--------------------------------------------------------------------------------
  CASE (16) !ppbar-> W+
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))
     FB=FB(5:-5:-1) !! inverse the quark order

    TMD_pair=&
    paramW_UD*(FA(2)*FB(-1)+FA(-1)*FB(2))&        !u*dbar+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(-3)*FB(2))&        !u*sbar+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(-5)*FB(2))&        !u*bbar+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(-1)*FB(4))&        !c*dbar+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(-3)*FB(4))&        !c*sbar+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(-5)*FB(4))        !c*bbar+bbar*c

!--------------------------------------------------------------------------------
  CASE (17) !ppbar-> W-
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))
     FB=FB(5:-5:-1) !! inverse the quark order

    TMD_pair=&
    paramW_UD*(FA(1)*FB(-2)+FA(-2)*FB(1))&        !d*ubar+ubar*d
    +paramW_US*(FA(3)*FB(-2)+FA(-2)*FB(3))&        !s*ubar+ubar*s
    +paramW_UB*(FA(5)*FB(-2)+FA(-2)*FB(5))&        !b*ubar+ubar*b
    +paramW_CD*(FA(1)*FB(-4)+FA(-4)*FB(1))&        !d*cbar+cbar*d
    +paramW_CS*(FA(3)*FB(-4)+FA(-4)*FB(3))&        !s*cbar+cbar*s
    +paramW_CB*(FA(5)*FB(-4)+FA(-4)*FB(5))        !b*cbar+cbar*b

!--------------------------------------------------------------------------------
  CASE (18) !ppbar-> W+ + W-
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(2))
     FB=FB(5:-5:-1) !! inverse the quark order

    TMD_pair=&
    paramW_UD*(FA(2)*FB(-1)+FA(1)*FB(-2)+FA(-2)*FB(1)+FA(-1)*FB(2))&    !u*dbar+d*ubar+ubar*d+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(3)*FB(-2)+FA(-2)*FB(3)+FA(-3)*FB(2))&    !u*sbar+s*ubar+ubar*s+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(5)*FB(-2)+FA(-2)*FB(5)+FA(-5)*FB(2))&    !u*bbar+b*ubar+ubar*b+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(1)*FB(-4)+FA(-4)*FB(1)+FA(-1)*FB(4))&    !c*dbar+d*cbar+cbar*d+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(3)*FB(-4)+FA(-4)*FB(3)+FA(-3)*FB(4))&    !c*sbar+s*cbar+cbar*s+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(5)*FB(-4)+FA(-4)*FB(5)+FA(-5)*FB(4))    !c*bbar+b*cbar+cbar*b+bbar*c
!--------------------------------------------------------------------------------
  CASE(20) !pp -> Higgs (unpol.part+lin.pol.part)
    FA=uTMDPDF_kT_50(x1,sqrt(k1),mu,Q2,process(2))
    FB=uTMDPDF_kT_50(x2,sqrt(k2),mu,Q2,process(2))
    TMD_pair=FA(0)*FB(0) !!!! unpolarized part

    FA=lpTMDPDF_kT_50(x1,sqrt(k1),mu,Q2,process(2))
    FB=lpTMDPDF_kT_50(x2,sqrt(k2),mu,Q2,process(2))
    TMD_pair=TMD_pair+FA(0)*FB(0) !!!! linearly polarized part
!--------------------------------------------------------------------------------
  CASE(21) !pp -> Higgs (unpol.part)
    FA=uTMDPDF_kT_50(x1,sqrt(k1),mu,Q2,process(2))
    FB=uTMDPDF_kT_50(x2,sqrt(k2),mu,Q2,process(2))
    TMD_pair=FA(0)*FB(0)

!--------------------------------------------------------------------------------
  CASE(22) !pp -> Higgs (lin.pol.part)
    FA=lpTMDPDF_kT_50(x1,sqrt(k1),mu,Q2,process(2))
    FB=lpTMDPDF_kT_50(x2,sqrt(k2),mu,Q2,process(2))
    TMD_pair=FA(0)*FB(0)

  !--------------------------------------------------------------------------------
  CASE (1001) !p+Cu->gamma* !!this is for E288
    FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
    FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(3))
    TMD_pair=116d0/567d0*(FA(2)*FB(-2)+FA(-2)*FB(2))+136d0/567d0*(FA(-2)*FB(1)+FA(2)*FB(-1))&
          +34d0/567d0*(FA(-1)*FB(2)+FA(1)*FB(-2))+29d0/567d0*(FA(-1)*FB(1)+FA(1)*FB(-1))&
          +1d0/9d0*(FA(-3)*FB(3)+FA(3)*FB(-3)+4d0*FA(-4)*FB(4)+4d0*FA(4)*FB(-4)+FA(-5)*FB(5)+FA(5)*FB(-5))

  !--------------------------------------------------------------------------------
  CASE (1002) !p+2H->gamma* !!this is for E772
    FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
    FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(3))
    TMD_pair=2d0/9d0*(FA(2)*FB(-2)+FA(-2)*FB(2))+2d0/9d0*(FA(-2)*FB(1)+FA(2)*FB(-1))&
          +1d0/18d0*(FA(-1)*FB(2)+FA(1)*FB(-2))+1d0/18d0*(FA(-1)*FB(1)+FA(1)*FB(-1))&
          +1d0/9d0*(FA(-3)*FB(3)+FA(3)*FB(-3)+4d0*FA(-4)*FB(4)+4d0*FA(4)*FB(-4)+FA(-5)*FB(5)+FA(5)*FB(-5))
  !--------------------------------------------------------------------------------
  CASE (1003) !pbar+W->gamma* !!this is for E537
    !Wolfram has A=183,    Z=74,    N=109
    FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
    FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(3))
    TMD_pair=296d0/1647d0*(FA(2)*FB(2)+FA(-2)*FB(-2))+436d0/1647d0*(FA(2)*FB(1)+FA(-2)*FB(-1))&
          +109d0/1647d0*(FA(1)*FB(2)+FA(-1)*FB(-2))+74d0/1647d0*(FA(1)*FB(1)+FA(-1)*FB(-1))&
          +1d0/9d0*(FA(3)*FB(3)+FA(-3)*FB(-3)+4d0*FA(4)*FB(4)+4d0*FA(-4)*FB(-4)+FA(5)*FB(5)+FA(-5)*FB(-5))
  !--------------------------------------------------------------------------------
  CASE (1004) !pminus+W->gamma* !!this is for E537
    !Wolfram has A=183,    Z=74,    N=109
    FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,process(2))
    FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,process(3))
    TMD_pair=296d0/1647d0*(FA(-2)*FB(2)+FA(2)*FB(-2))+436d0/1647d0*(FA(-2)*FB(1)+FA(2)*FB(-1))&
          +109d0/1647d0*(FA(-1)*FB(2)+FA(1)*FB(-2))+74d0/1647d0*(FA(-1)*FB(1)+FA(1)*FB(-1))&
          +1d0/9d0*(FA(-3)*FB(3)+FA(3)*FB(-3)+4d0*FA(-4)*FB(4)+4d0*FA(4)*FB(-4)+FA(-5)*FB(5)+FA(5)*FB(-5))

  CASE DEFAULT
    write(*,*) ErrorString('undefined process: ',moduleName),process
    write(*,*) color('Evaluation stop',c_red_bold)
    stop
 END SELECT


end function TMD_pair


!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!!! The hadron tensonr for the DY icludes Z + gamma, evaluated at FA and FB
function XTMD_pairForDYwithZgamma(FAB,Q2)
     real(dp)::XTMD_pairForDYwithZgamma,Q2
    !!cross-seciton parameters
     real(dp),dimension(-5:5):: FAB

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

     XTMD_pairForDYwithZgamma=&
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

end function XTMD_pairForDYwithZgamma
