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
  CASE (2,3,4,5,20,21,22) !Delta^{GG'}z_{+l}z_{+f}f1f1
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,h1)
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,-h2)!!! -h2, to multiply quarks by anti-quarks in FAB
     FAB=FA*FB

    TMD_pair=XTMD_pairZpZp(FAB,Q2)

  !--------------------------------------------------------------------------------
  CASE (6,7,23,24) !Delta^{GG'}z_{-l}z_{-f}{f1f1}_A
     FA=uTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,h1)
     FB=uTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,-h2)!!! -h2, to multiply quarks by anti-quarks in FAB
     FAB=FA*FB

    TMD_pair=XTMD_pairZmZm_A(FAB,Q2)

  !--------------------------------------------------------------------------------
  CASE (8,9,10,30,31,32) !Delta^{GG'}z_{+l}r_{+f}h1h1
     FA=BoerMuldersTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,h1)
     FB=BoerMuldersTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,-h2)!!! -h2, to myltiply quarks by anti-quarks in FAB
     FAB=FA*FB

    TMD_pair=XTMD_pairZpRp(FAB,Q2)

  !--------------------------------------------------------------------------------
  CASE (11,12,35,36) !Delta^{GG'}z_{+l}r_{-f}{h1h1}_A
     FA=BoerMuldersTMDPDF_kT_5(x1,sqrt(k1),mu,Q2,h1)
     FB=BoerMuldersTMDPDF_kT_5(x2,sqrt(k2),mu,Q2,-h2)!!! -h2, to myltiply quarks by anti-quarks in FAB
     FAB=FA*FB

    TMD_pair=XTMD_pairZpRm_A(FAB,Q2)


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

     XTMD_pairZpZp=&
     zP_gg_L*(&!gamma-part
      zP_gg_U*FAB(2)&
      +zP_gg_D*FAB(1)&
      +zP_gg_S*FAB(3)&
      +zP_gg_C*FAB(4)&
      +zP_gg_B*FAB(5)&
      +zP_gg_U*FAB(-2)&
      +zP_gg_D*FAB(-1)&
      +zP_gg_S*FAB(-3)&
      +zP_gg_C*FAB(-4)&
      +zP_gg_B*FAB(-5))&
     +&!gamma-Z interference
     zP_gZ_L*(&
      zP_gZ_U*FAB(2)&
      +zP_gZ_D*FAB(1)&
      +zP_gZ_S*FAB(3)&
      +zP_gZ_C*FAB(4)&
      +zP_gZ_B*FAB(5)&
      +zP_gZ_U*FAB(-2)&
      +zP_gZ_D*FAB(-1)&
      +zP_gZ_S*FAB(-3)&
      +zP_gZ_C*FAB(-4)&
      +zP_gZ_B*FAB(-5))*&
      2d0*Q2*(Q2-MZ2)/((Q2-MZ2)**2+GammaZ2*MZ2)&
     +&!ZZ-contributions
     zP_ZZ_L*(&
       zP_ZZ_U*FAB(2)&
      +zP_ZZ_D*FAB(1)&
      +zP_ZZ_S*FAB(3)&
      +zP_ZZ_C*FAB(4)&
      +zP_ZZ_B*FAB(5)&
      +zP_ZZ_U*FAB(-2)&
      +zP_ZZ_D*FAB(-1)&
      +zP_ZZ_S*FAB(-3)&
      +zP_ZZ_C*FAB(-4)&
      +zP_ZZ_B*FAB(-5))*&
      Q2*Q2/((Q2-MZ2)**2+GammaZ2*MZ2)

end function XTMD_pairZpZp

!!! Combination Delta^{GG'} z_{+l}r_{+f} FF
function XTMD_pairZpRp(FAB,Q2)
     real(dp)::XTMD_pairZpRp
     real(dp),intent(in)::Q2
    !!cross-seciton parameters
     real(dp),dimension(-5:5),intent(in):: FAB


     XTMD_pairZpRp=&
     zP_gg_L*(&!gamma-part   !!!!rP_gg=zP_gg
      zP_gg_U*FAB(2)&
      +zP_gg_D*FAB(1)&
      +zP_gg_S*FAB(3)&
      +zP_gg_C*FAB(4)&
      +zP_gg_B*FAB(5)&
      +zP_gg_U*FAB(-2)&
      +zP_gg_D*FAB(-1)&
      +zP_gg_S*FAB(-3)&
      +zP_gg_C*FAB(-4)&
      +zP_gg_B*FAB(-5))&
     +&!gamma-Z interference  !!!!rP_gZ=zP_gZ
     zP_gZ_L*(&
      zP_gZ_U*FAB(2)&
      +zP_gZ_D*FAB(1)&
      +zP_gZ_S*FAB(3)&
      +zP_gZ_C*FAB(4)&
      +zP_gZ_B*FAB(5)&
      +zP_gZ_U*FAB(-2)&
      +zP_gZ_D*FAB(-1)&
      +zP_gZ_S*FAB(-3)&
      +zP_gZ_C*FAB(-4)&
      +zP_gZ_B*FAB(-5))*&
      2d0*Q2*(Q2-MZ2)/((Q2-MZ2)**2+GammaZ2*MZ2)&
     +&!ZZ-contributions
     zP_ZZ_L*(&
       rP_ZZ_U*FAB(2)&
      +rP_ZZ_D*FAB(1)&
      +rP_ZZ_S*FAB(3)&
      +rP_ZZ_C*FAB(4)&
      +rP_ZZ_B*FAB(5)&
      +rP_ZZ_U*FAB(-2)&
      +rP_ZZ_D*FAB(-1)&
      +rP_ZZ_S*FAB(-3)&
      +rP_ZZ_C*FAB(-4)&
      +rP_ZZ_B*FAB(-5))*&
      Q2*Q2/((Q2-MZ2)**2+GammaZ2*MZ2)

end function XTMD_pairZpRp

!!! Combination Delta^{GG'} z_{-l}z_{-f} {FF}_A
function XTMD_pairZmZm_A(FAB,Q2)
     real(dp)::XTMD_pairZmZm_A
     real(dp),intent(in)::Q2
    !!cross-seciton parameters
     real(dp),dimension(-5:5),intent(in):: FAB

     XTMD_pairZmZm_A=&  !zM_gg =0
     zM_gZ_L*(& !gamma-Z interference
      zM_gZ_U*FAB(2)&
      +zM_gZ_D*FAB(1)&
      +zM_gZ_S*FAB(3)&
      +zM_gZ_C*FAB(4)&
      +zM_gZ_B*FAB(5)&
      -zM_gZ_U*FAB(-2)&
      -zM_gZ_D*FAB(-1)&
      -zM_gZ_S*FAB(-3)&
      -zM_gZ_C*FAB(-4)&
      -zM_gZ_B*FAB(-5))*&
      2d0*Q2*(Q2-MZ2)/((Q2-MZ2)**2+GammaZ2*MZ2)&
     +&!ZZ-contributions
       zM_ZZ_L*(&
       zM_ZZ_U*FAB(2)&
      +zM_ZZ_D*FAB(1)&
      +zM_ZZ_S*FAB(3)&
      +zM_ZZ_C*FAB(4)&
      +zM_ZZ_B*FAB(5)&
      -zM_ZZ_U*FAB(-2)&
      -zM_ZZ_D*FAB(-1)&
      -zM_ZZ_S*FAB(-3)&
      -zM_ZZ_C*FAB(-4)&
      -zM_ZZ_B*FAB(-5))*&
      Q2*Q2/((Q2-MZ2)**2+GammaZ2*MZ2)

end function XTMD_pairZmZm_A

!!! Combination Delta^{GG'} i z_{+l}r_{-f} {FF}_A
function XTMD_pairZpRm_A(FAB,Q2)
     real(dp)::XTMD_pairZpRm_A
     real(dp),intent(in)::Q2
    !!cross-seciton parameters
     real(dp),dimension(-5:5),intent(in):: FAB

     XTMD_pairZpRm_A=&  !zM_gg =0
     zP_gZ_L*(& !gamma-Z interference
      rM_gZ_U*FAB(2)&
      +rM_gZ_D*FAB(1)&
      +rM_gZ_S*FAB(3)&
      +rM_gZ_C*FAB(4)&
      +rM_gZ_B*FAB(5)&
      -rM_gZ_U*FAB(-2)&
      -rM_gZ_D*FAB(-1)&
      -rM_gZ_S*FAB(-3)&
      -rM_gZ_C*FAB(-4)&
      -rM_gZ_B*FAB(-5))*&
      2d0*Q2*sqrt(MZ2*GammaZ2)/((Q2-MZ2)**2+GammaZ2*MZ2)
      !! no ZZ-term


end function XTMD_pairZpRm_A
