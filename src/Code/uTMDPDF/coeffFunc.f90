!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of uTMDPDF module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for the mathing coefficient--------------------------------------
!!!--------------The order is accumulative pertrubative order of coefficient =0,1,2 (LO,NLO,NNLO)---------
!!!-------------------------------------------------------------------------------------------------------

  !!! the function which contains the functions of parameterizations
  function parametrizationString(z)
  real(dp)::z,lz,l1z,zzz
  real(dp),dimension(1:parametrizationLength)::parametrizationString

  real(dp),dimension(0:23)::t !!! chebhyshev polynomial
  integer::i

      zzz=2*z-1
      t(0)=1_dp
      t(1)=zzz
      do i=2,23
        t(i)=2*zzz*t(i-1)-t(i-2)
      end do

      lz=Log(z)
      l1z=Log(1d0-z)

      parametrizationString=(/&
        l1z,l1z**2,l1z**3,&
	    1d0/z,lz/z,&
	    lz,lz**2,lz**3,&
	    1d0,z,z**2,&
	    z/(1d0-z)*lz, z*lz, (z**2)*lz,&
	    z/(1d0-z)*lz**2,z*lz**2,&
	    (lz/(1d0-z)+1d0)*l1z,lz*l1z,z*lz*l1z,&
	    (1d0-z)/z*l1z, (1d0-z)*l1z, ((1d0-z)**2)*l1z,(1d0-z)*l1z**2/)
  
  end function parametrizationString
  
  !!! the function which contains 
    !!! int_z^1 parameterization at values of z -> 1
    !!! it is used to estimate integration error at z~1
  function parametrizationStringAt1(z)
  real(dp)::z,l1z
  real(dp),dimension(1:parametrizationLength)::parametrizationStringAt1
    l1z=Log(1d0-z)
     parametrizationStringAt1=(/(1d0-z)*(l1z-1),(1d0-z)*(2d0-2d0*l1z+l1z**2),&
	      (1d0-z)*(-6d0+6d0*l1z-3d0*l1z**2+l1z**3),&
	      0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,&
	      ((1d0-z)**2)*l1z,&
	      0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  
  end function parametrizationStringAt1

  
  !!!!Each coefficient is split to delta, sing x->1, regular
    
  !!!!!coefficient function q<-q delta-part
pure function C_q_q_delta(alpha,Nf,Lmu)
  real(dp),intent(in)::alpha,Nf,Lmu
  real(dp)::C_q_q_delta
  
 C_q_q_delta=1d0
 
  if(order_global>=1) then
      !!!checked (3.12.17) AV!!!
      C_q_q_delta=C_q_q_delta+alpha*(-4d0/3d0*zeta2-4d0*Lmu)

  if(order_global>=2) then
      !!!checked (3.12.17) AV!!!
     C_q_q_delta=C_q_q_delta+alpha*alpha*(&
     -2416d0/81d0-134d0/3d0*zeta2+448d0/9d0*zeta3+200d0/9d0*zeta4+&
     Nf*(352d0/243d0+20d0/9d0*zeta2+56d0/27d0*zeta3)+&
     Lmu*(-14d0-140d0/3d0*zeta2+16d0/3d0*zeta3+Nf*(4d0/9d0+40d0/9d0*zeta2))+&
     Lmu**2*(-14d0+4d0/3d0*Nf-128d0/9d0*zeta2))

  if(order_global>=3) then
    C_q_q_delta=C_q_q_delta+alpha**3*(&
    Nf**2*(2800d0/19683d0 - (496d0*zeta2)/81d0 - (3712d0*zeta3)/729d0 - (88d0*zeta4)/81d0)&
    +Nf*(212644d0/6561d0 + (224116d0*zeta2)/729d0 + (83452d0*zeta3)/729d0 - (5960d0*zeta2*zeta3)/81d0 &
      - (1988d0*zeta4)/81d0 - (8144d0*zeta5)/81d0)&
    -578966d0/2187d0 - (667234d0*zeta2)/243d0 + (13954d0*zeta3)/81d0 + (30212d0*zeta2*zeta3)/27d0 &
    -(244d0*zeta3**2)/3d0 + (12796d0*zeta4)/9d0  - (1576d0*zeta5)/3d0 - (59468d0*zeta6)/81d0&
    +Lmu*(&
        Nf**2*(428d0/729d0 - (560d0*zeta2)/81d0 - (32d0*zeta3)/81d0)&
        + Nf*(6496d0/243d0 + (80536d0*zeta2)/243d0 - (112d0*zeta3)/3d0 - (1064d0*zeta4)/9d0)&
        -5947d0/27d0 - (232612d0*zeta2)/81d0 + (35200d0*zeta3)/27d0 + (15680d0*zeta2*zeta3)/27d0&
        + (3364d0*zeta4)/3d0 - (4960d0*zeta5)/9d0)&
    +Lmu**2*(-302d0 + Nf**2*(-8d0/27d0 - (80d0*zeta2)/27d0) - (26476d0*zeta2)/27d0 &
        + Nf*(340d0/9d0 + (7744d0*zeta2)/81d0 - (32d0*zeta3)/9d0) + (112d0*zeta3)/3d0 + (12800d0*zeta4)/27d0)&
    +Lmu**3*(-84d0 - 16d0*Nf**2/27d0 - 896d0*zeta2/9d0 + Nf*(128d0/9d0 + 256d0*zeta2/27d0) - 4096d0*zeta3/81d0)&
    )

  end if
  end if
  end if

end function C_q_q_delta
  
  !!!!!coefficient function g<-g delta-part
pure function C_g_g_delta(alpha,Nf,Lmu)
  real(dp),intent(in)::alpha,Nf,Lmu
  real(dp)::C_g_g_delta
  
  C_g_g_delta=1d0
 
  if(order_global>=1) then
  !!!checked (20.11.17) AV!!!
      C_g_g_delta=C_g_g_delta+alpha*(-3d0*zeta2-(11d0-2d0/3d0*Nf)*Lmu)

  if(order_global>=2) then
  !!!checked (23.06.22) AV!!!
     C_g_g_delta=C_g_g_delta+alpha*alpha*(&
    -112d0 - 56d0*Nf**2/81d0 - 201d0*zeta2/2d0 - 72d0*Lmu**2*zeta2 +&
    Lmu*(-96d0 + 32d0*Nf/3d0 - 108d0*zeta3) + Nf*(548d0/27 + 5d0*zeta2 - 28d0*zeta3/3d0) + 154d0*zeta3 + 225d0*zeta4/4d0)

  if(order_global>=3) then
  C_g_g_delta=C_g_g_delta+alpha**3*(&
    Nf**3*(-752d0/2187d0 + (16d0*zeta3)/27d0)&
    + Nf**2*(-73577d0/2187d0 - (200d0*zeta2)/81d0 - (368d0*zeta3)/81d0 - (20d0*zeta4)/9d0)&
    + Nf*(1033259d0/1458d0 + (27305d0*zeta2)/81d0 - (17762d0*zeta3)/81d0 + (122d0*zeta2*zeta3)/3d0 &
    - (2767d0*zeta4)/18d0 + (556d0*zeta5)/9d0)&
    -698456d0/243d0 - (213865d0*zeta2)/54d0 + (1489d0*zeta3)/9d0 + 429d0*zeta2*zeta3 + 2337d0*zeta3**2&
    + (15395d0*zeta4)/4d0 - 2046d0*zeta5  - (28783d0*zeta6)/16d0&
    +Lmu*((112d0*Nf**3)/243d0&
      +Nf**2*(-4471d0/162d0 - (10d0*zeta2)/3d0 + (56d0*zeta3)/9d0)&
      +Nf*(25175d0/54d0 + (904d0*zeta2)/3d0 + (104d0*zeta3)/3d0 - (15d0*zeta4)/2d0)&
      -4597d0/2d0 - (8855d0*zeta2)/2d0 - 3130d0*zeta3 + 3780d0*zeta2*zeta3 + (495d0*zeta4)/4d0 + 2160d0*zeta5)&
    + Lmu**2*(-561d0 - (38d0*Nf**2)/9d0 - 3216d0*zeta2 + Nf*(311d0/3d0 + 160d0*zeta2) + 2700d0*zeta4)&
    - Lmu**3*576d0*zeta3&
    )

  end if
  end if
  end if
  
end function C_g_g_delta
  
  !!!!!coefficient function q<-q singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
subroutine Set_CoeffSing1_q_q(alpha,Nf,Lmu)
  real(dp),intent(in)::Nf,alpha,Lmu
  real(dp)::s1,s2,s3
    
  s1=0d0!!!coeff 1/(1-x)
  s2=0d0!!!coeff log(1-x)/(1-x)
  s3=0d0!!!coeff log(1-x)^2/(1-x)
  

  if(order_global>=1) then    
    !!!checked (20.11.17) AV!!!
    s1=s1+alpha*(-16d0/3d0)*Lmu

  if(order_global>=2) then
    !!!checked (21.11.17) AV!!!
    s1=s1+alpha*alpha*(&
       -3232d0/27d0+112d0*zeta3+448d0/81d0*Nf+&
       Lmu*(-1072d0/9d0+352d0/9d0*zeta2+160d0/27d0*Nf)+&
       Lmu**2*(-8d0+16d0/9d0*Nf))

     s2=s2+alpha*alpha*(256d0/9d0*(Lmu**2))

  if(order_global>=3) then
    s1=s1+alpha**3*(&
      Nf**2*(-7424d0/2187d0 - (128d0*zeta3)/27d0)&
      + Nf*(332632d0/729d0 - (11680d0*zeta2)/243d0 - (15712d0*zeta3)/81d0 - (16d0*zeta4)/9d0)&
      -1188116d0/243d0 + (89632d0*zeta2)/81d0 + (49312d0*zeta3)/9d0 - (2560d0*zeta2*zeta3)/3d0&
      + 616d0*zeta4 - 2304d0*zeta5&
      +Lmu*(&
        - (1600d0*Nf**2)/243d0&
        + Nf*(321104d0/729d0 - (7360d0*zeta2)/81d0 - (5504d0*zeta3)/81d0)&
        -961208d0/243d0  + (49312d0*zeta2)/27d0 + (37760d0*zeta3)/27d0 - (34592d0*zeta4)/27d0)&
      +Lmu**2*(&
        -9280d0/9d0 - (320d0*Nf**2)/81d0 + Nf*(4112d0/27d0 - (1280d0*zeta2)/27d0) + 512d0*zeta2 - (256d0*zeta3)/9d0)&
      +Lmu**3*(&
        -208d0/9d0 + (320d0*Nf)/27d0 - (64d0*Nf**2)/81d0 + (2048d0*zeta2)/27d0)&
      )

    s2=s2+alpha**3*(&
      Lmu**3*(1792d0/9d0 - (512d0*Nf)/27d0) &
      + Lmu**2*(34304d0/27d0 - (5120d0*Nf)/81d0 - (10240d0*zeta2)/27d0) &
      + Lmu*(103424d0/81d0 - (14336d0*Nf)/243d0 - (3584d0*zeta3)/3d0)&
      )

    s3=s3+alpha**3*(-2048d0/27d0*Lmu**3)

  end if
  end if
  end if

  CoeffSing1_q_q=(/s1,s2,s3/)
  
end subroutine Set_CoeffSing1_q_q
  
  !!!!!coefficient function g<-g singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
subroutine Set_CoeffSing1_g_g(alpha,Nf,Lmu)
  real(dp),intent(in)::Nf,alpha,Lmu
  real(dp)::s1,s2,s3
    
  s1=0d0!!!coeff 1/(1-x)
  s2=0d0!!!coeff log(1-x)/(1-x)
  s3=0d0!!!coeff log(1-x)^2/(1-x)
 
  if(order_global>=1) then  
  !!!checked (20.11.17) AV!!!
    s1=s1+alpha*(-12d0)*Lmu

  if(order_global>=2) then
  !!!checked (23.06.22) AV!!!
    s1=s1+alpha*alpha*(&
    -808d0/3d0 + Lmu**2*(66d0 - 4d0*Nf) + 112d0*Nf/9d0 + Lmu*(-268d0 + 40d0*Nf/3d0 + 108d0*zeta2) + 252d0*zeta3)
    
     s2=s2+alpha*alpha*144d0*(Lmu**2)

  if(order_global>=3) then
    s1=s1+alpha**3*(&
    Nf**2*(-1856d0/243d0 - (32d0*zeta3)/3d0) &
    +Nf*(83158d0/81d0 - (1160d0*zeta2)/9d0 - (3928d0*zeta3)/9d0 - 4d0*zeta4) &
    -297029d0/27d0 + (8816d0*zeta2)/3d0 + 12328d0*zeta3 - 2340d0*zeta2*zeta3 + 1386d0*zeta4 - 5184d0*zeta5 &
    +Lmu*(-18086d0/3d0 + (16d0*Nf**2)/9d0 + 5226d0*zeta2 + 132d0*zeta3 + &
      Nf*(4484d0/9d0 - 260d0*zeta2 + 152d0*zeta3) - 3591d0*zeta4) &
    +Lmu**2*(540d0 + Nf*(-52d0 - 12d0*zeta2) + 198d0*zeta2 + 1296d0*zeta3) &
    +Lmu**3*(242d0 - (88d0*Nf)/3d0 + (8d0*Nf**2)/9d0 + 864d0*zeta2) &
    )

     s2=s2+alpha**3*(&
     Lmu**2*(6432d0 - 320d0*Nf - 2160d0*zeta2) + Lmu*(6464d0 - (896d0*Nf)/3 - 6048d0*zeta3))

     s3=s3+alpha**3*(-864d0*Lmu**3)
  end if
  end if
  end if
  
  CoeffSing1_g_g=(/s1,s2,s3/)
  
end subroutine Set_CoeffSing1_g_g
  
  !!!!!coefficient function q<-q regular-part  
  subroutine Set_Coeff_q_q(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_q=(/&
  0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
  0d0,0d0,&     !1/x, Log(x)/x
  0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
  0d0,0d0,0d0,& !1 ,x , x^2
  0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
  0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
  0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
  0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  if(order_global>=1) then
          
      Coeff_q_q=Coeff_q_q+alpha*(/&
	0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
	0d0,0d0,&     !1/x, Log(x)/x
	0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
	8d0/3d0*(1d0+Lmu),8d0/3d0*(Lmu-1d0),0d0,& !1 ,x , x^2
	0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
	0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
	0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
	0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
    
     Coeff_q_q=Coeff_q_q+alpha*alpha*((/&!!!exact part
      200d0/9d0 - 256d0*Lmu/9d0 - 256d0*Lmu**2/9d0, -64d0/9d0, 0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      0d0,0d0,&     !1/x, Log(x)/x
      -8d0 - 32d0*Lmu**2/9d0 + 40d0*Nf/27d0 + Lmu*(-152d0/9d0 + 16d0*Nf/9d0), &
	-2d0 - 40d0*Lmu/9d0 + 4d0*Nf/9d0, -20d0/27d0,& !Log(x),Log^2(x),Log^3(x)
      Lmu**2*(-28d0/9d0 - 8d0*Nf/9d0) - 296d0*Nf/81d0 + Lmu*(64d0/3d0 - 32d0*Nf/27d0 - 176d0*zeta2/9d0), &
	Lmu**2*(100d0/9d0 - 8d0*Nf/9d0) - 152d0*Nf/81d0 + Lmu*(880d0/9d0 - 128d0*Nf/27d0 - 176d0*zeta2/9d0), 0d0,& !1 ,x , x^2
      -128d0*Lmu**2/9d0 + 80d0*Nf/27d0 + Lmu*(-112d0/3d0 + 32d0*Nf/9d0), &
	32d0*Lmu**2/3d0 + Lmu*(184d0/9d0 - 16d0*Nf/9d0) - 40d0*Nf/27d0, 0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      -16d0*Lmu + 8d0*Nf/9d0, 104d0*Lmu/9d0 - 4d0*Nf/9d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      256d0*Lmu/9d0, -128d0*Lmu/9d0, -128d0*Lmu/9d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      0d0, -128d0*Lmu/9d0 + 128d0*Lmu**2/9d0, 0d0, 0d0/)& !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
    +(/&!!aprroximate 
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      0d0,0d0,&     !1/x, Log(x)/x
      0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
      1076.6744297226016d0, 7792.719665777814d0, 111.49810429898287d0,& !1 ,x , x^2
      8980.334190376141d0, -3795.008745809993d0, 82.30795871692112d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      -201.0129463471822d0, 206.75145891009598d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      5603.371344939401d0, -526.1352578350587d0, 1382.8610999663256d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      1092.9256332669593d0, 2547.784733022028d0, -147.17479558391307d0, 3.564983084988843d0/)) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
    end if
    
!  write(*,*) 'regularPart=', regularPart/x
  end if
  end subroutine Set_Coeff_q_q
 
  !!!!!coefficient function q<-g regular-part  
  subroutine Set_Coeff_q_g(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_g=(/&
  0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
  0d0,0d0,&     !1/x, Log(x)/x
  0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
  0d0,0d0,0d0,& !1 ,x , x^2
  0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
  0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
  0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
  0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  if(order_global>=1) then
    Coeff_q_g=Coeff_q_g+alpha*(/ &
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      0d0,0d0,&     !1/x, Log(x)/x
      0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
      -Lmu,2d0+2d0*Lmu,-2d0-2d0*Lmu,& !1 ,x , x^2
      0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  
    if(order_global>=2) then
     Coeff_q_g=Coeff_q_g+alpha*alpha*((/&!!!Exact part
      -5d0/3d0 + 26d0*Lmu**2/3d0, 10d0*Lmu/3d0, 5d0/9d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      172d0/9d0 - 52d0*Lmu/3d0 + 4d0*Lmu**2 - 8d0*zeta2, 0d0,&     !1/x, Log(x)/x
      34d0/3d0 - 22d0*Lmu/3d0 + 14d0*Lmu**2/3d0, -7d0/6d0 + 14d0*Lmu/3d0, 7d0/9d0,& !Log(x),Log^2(x),Log^3(x)
      7d0*Lmu**2/3d0, 80d0*Lmu**2/3d0, -31d0*Lmu**2,& !1 ,x , x^2
      0d0, 80d0*Lmu**2/3d0, -16d0*Lmu**2/3d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      0d0,-52d0*Lmu**2/3d0, 52d0*Lmu**2/3d0,0d0/)& !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
      +(/&!!!Approximate part
	0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
	0d0,0d0,&     !1/x, Log(x)/x
	0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
	11804.917158263002d0 - 241.5977366488703d0*Lmu, &
	  43420.91122559322d0 - 1282.2433019065484d0*Lmu, &
	  -3972.3493167408856d0 - 200.38789519756068d0*Lmu,& !1 ,x , x^2
	51274.90217179903d0 - 1729.992223130878d0*Lmu, &
	  -11322.606629412014d0 + 831.0967099095751d0*Lmu, &
	  1490.8267572457225d0 - 43.62058264210284d0*Lmu,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
	365.38010578193774d0 + 17.18487164919196d0*Lmu, &
	  -402.6162297779727d0 - 2.712778270327616d0*Lmu,&     !x/(1-x)Log^2(x), x Log^2(x)
	21027.101170409886d0 - 1172.7172087577621d0*Lmu, &
	  19269.1080641264d0 + 259.31259625889305d0*Lmu, &
	  11357.814388977074d0 - 185.71439659057268d0*Lmu,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
	11798.406220601846d0 - 263.2306304278585d0*Lmu, &
	  29476.38189858057d0 - 264.4846500672005d0*Lmu, &
	  -227.4273932698467d0 + 82.74615814331906d0*Lmu, &
	  15.408922558612357d0 - 6.5920606591307065d0*Lmu/)) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
    end if
    
  end if
  end subroutine Set_Coeff_q_g
  
    !!!!!coefficient function g<-q regular-part  
  subroutine Set_Coeff_g_q(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_q=(/&
    0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
    0d0,0d0,&     !1/x, Log(x)/x
    0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
    0d0,0d0,0d0,& !1 ,x , x^2
    0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
    0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
    0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
    0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  if(order_global>=1) then
    Coeff_g_q=Coeff_g_q+alpha*(/&
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      -16d0/3d0*Lmu,0d0,&     !1/x, Log(x)/x
      0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
      16d0/3d0*Lmu,8d0/3d0*(1d0-Lmu),0d0,& !1 ,x , x^2
      0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
     Coeff_g_q=Coeff_g_q+alpha*alpha*((/ &!!Exact part
      -184d0/9d0 + 208d0*Lmu**2/9d0 + 32d0*Nf/27d0 + Lmu*(-208d0/3d0 + 32d0*Nf/9d0), &
	 -44d0/9d0 - 80d0*Lmu/9d0 + 8d0*Nf/9d0, -40d0/27d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      -12640d0/27d0 - 248d0*Lmu**2/3d0 + 896d0*Nf/81d0 + 352d0*zeta2/3d0 + Lmu*(-16d0 + 320d0*Nf/27d0 + 16d0*zeta2)&
	  + 192d0*zeta3, -32d0*Lmu**2,&     !1/x, Log(x)/x
      -200d0/3d0 + 448d0*Lmu/9d0 - 224d0*Lmu**2/9d0, 112d0/9d0 - 224d0*Lmu/9d0, -112d0/27d0,& !Log(x),Log^2(x),Log^3(x)
      640d0*Lmu**2/9d0 - 320d0*Lmu*Nf/27d0, 56d0*Lmu**2/9d0 + 208d0*Lmu*Nf/27d0, 32d0*Lmu**2/3d0,& !1 ,x , x^2
      0d0,-320d0*Lmu**2/9d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      416d0*Lmu**2/9d0 + 64d0*Lmu*Nf/9d0, -208d0*Lmu**2/9d0 - 32d0*Lmu*Nf/9d0,0d0,0d0/)&
      +(/&!!approximate part
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      0d0,0d0,&     !1/x, Log(x)/x
      0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
      25387.766684267783d0 + 995.0373392459464d0*Lmu - 155.5078210719438d0*Nf, &
	128641.80191589122d0 + 3632.526455640396d0*Lmu - 529.7740230354465d0*Nf, &
	-4691.90206636766d0 - 626.2634197396076d0*Lmu + 33.2493657381115d0*Nf,& !1 ,x , x^2
      149304.47081486747d0 + 4104.585347056393d0*Lmu - 643.5386512087849d0*Nf, &
	-49319.02516115282d0 - 501.0681979811703d0*Lmu + 160.4570922525814d0*Nf, &
	3077.798040737399d0 + 226.19867367258644d0*Lmu - 15.139140368040819d0*Nf,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      -1297.381684493681d0 - 5.587857511419472d0*Lmu + 2.266359003377411d0*Nf, &
	1379.187692679483d0 - 13.897150833036367d0*Lmu - 2.2865882004537834d0*Nf,&     !x/(1-x)Log^2(x), x Log^2(x)
      79157.51026314059d0 + 1713.027825100937d0*Lmu - 313.45913506070985d0*Nf, &
	20294.80059647656d0 + 1827.284781273226d0*Lmu - 169.33020796238313d0*Nf, &
	27727.02916158908d0 + 870.6474261030947d0*Lmu - 115.18065110064754d0*Nf,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      25431.533693585556d0 + 1011.076509038249d0*Lmu - 138.52016675095615d0*Nf, &
	62304.61943144707d0 + 2503.327007358015d0*Lmu - 297.70931938249043d0*Nf,&
	-1802.796656932475d0 - 3.5567918590974164d0*Lmu + 4.76862827955143d0*Nf, &
	5.352770139788647d0 - 9.08059951934387d0*Lmu + 0.9080131334622416d0*Nf/))
    end if
    
  end if
  end subroutine Set_Coeff_g_q
  
      !!!!!coefficient function g<-g regular-part  
  subroutine Set_Coeff_g_g(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_g=(/&
    0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
    0d0,0d0,&     !1/x, Log(x)/x
    0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
    0d0,0d0,0d0,& !1 ,x , x^2
    0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
    0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
    0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
    0d0,0d0,0d0,0d0/)
  if(order_global>=1) then
    Coeff_g_g=Coeff_g_g+alpha*Lmu*(/&
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      -12d0,0d0,&     !1/x, Log(x)/x
      0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
      24d0,-12d0,12d0,& !1 ,x , x^2
      0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      0d0,0d0,0d0,0d0/)
  
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
     Coeff_g_g=Coeff_g_g+alpha*alpha*((/ &!!Exact part
      6d0 - 144d0*Lmu - 144d0*Lmu**2 - 2d0*Nf, -36d0, 0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      -3160d0/3d0 + Lmu**2*(-198d0 - 4d0*Nf/9d0) + 226d0*Nf/9d0 + 264d0*zeta2 + Lmu*(244d0*Nf/9d0 + 36d0*zeta2)&
	+ 432d0*zeta3,-72d0*Lmu**2,&     !1/x, Log(x)/x
      -293d0 + 74d0*Nf/3d0 + Lmu**2*(-72d0 + 16d0*Nf/3d0) + Lmu*(12d0 + 24d0*Nf), &
	3d0 + 6d0*Nf + Lmu*(-72d0 + 16d0*Nf/3d0), -12d0 + 8d0*Nf/9d0,& !Log(x),Log^2(x),Log^3(x)
      -4d0*Lmu*Nf/3d0 + Lmu**2*(84d0 + 32d0*Nf/3d0), Lmu**2*(-150d0 - 20d0*Nf/3d0) - 4d0*Lmu*Nf/3d0, &
	Lmu**2*(198d0 + 4d0*Nf/9d0) - 340d0*Lmu*Nf/9d0,& !1 ,x , x^2
      -72d0*Lmu**2, 24d0*Lmu*Nf + Lmu**2*(-216d0 + 16d0*Nf/3d0), 72d0*Lmu**2,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      0d0,16d0*Lmu*Nf/3d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      144d0*Lmu**2, 144d0*Lmu**2, -144d0*Lmu**2,0d0/)&
      +(/&!!approximate part
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      0d0,0d0,&     !1/x, Log(x)/x
      0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
      28350.905925013587d0 - 11826.30667028584d0*Lmu - 1165.597999487027d0*Nf, &
	204106.38400427837d0 - 83819.2024142788d0*Lmu - 9112.918558788755d0*Nf, &
	-3654.296204682914d0 - 4654.703280638828d0*Lmu - 202.82515955533333d0*Nf,& !1 ,x , x^2
      228606.3613205384d0 - 100331.34185957731d0*Lmu - 10439.786162275559d0*Nf, &
	-92515.7827039022d0+ 44024.02341566855d0*Lmu + 4572.107887712472d0*Nf, &
	4598.901824247715d0 + 129.71334826178415d0*Lmu - 82.57482307730146d0*Nf,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      -5961.667327556318d0 + 539.0169565557716d0*Lmu + 242.87579888233213d0*Nf, &
	6430.443013955799d0 - 686.8969602769131d0*Lmu - 258.29777371318767d0*Nf,&     !x/(1-x)Log^2(x), x Log^2(x)
      144229.2381370356d0 - 56326.8843090263d0*Lmu - 6707.548773597168d0*Nf, &
	-11622.452033932963d0 + 209.23681132222114d0*Lmu + 958.7180437906424d0*Nf, &
	35791.91490687339d0 - 18745.132737698765d0*Lmu - 1542.1541360841081d0*Nf,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      28788.30003557479d0 - 11546.565659033551d0*Lmu - 1171.3757772648048d0*Nf, &
	67492.96763028382d0 - 35083.678877258644d0*Lmu - 2764.3635575309518d0*Nf, &
	-3881.2814494298837d0 + 1923.6823271680485d0*Lmu + 176.34274530829245d0*Nf, &
	18.352295609756112d0 - 0.22756238102609413d0*Lmu - 0.04234942977167983d0*Nf/))
    end if
    
  end if
  end subroutine Set_Coeff_g_g

   !!!!!coefficient function q<-qb regular-part  
  subroutine Set_Coeff_q_qb(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qb=(/&
    0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
    0d0,0d0,&     !1/x, Log(x)/x
    0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
    0d0,0d0,0d0,& !1 ,x , x^2
    0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
    0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
    0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
    0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
     Coeff_q_qb=Coeff_q_qb+alpha*alpha*(/&
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      0d0,0d0,&     !1/x, Log(x)/x
      -4d0/3d0 + 16d0*Lmu/9d0, 8d0*Lmu/9d0, 4d0/27d0,& !Log(x),Log^2(x),Log^3(x)
      540.6779037661292d0 + 488.677462854424d0*Lmu, &
	4561.881499265996d0 + 3844.6673748808585d0*Lmu, &
	330.3186826846845d0 + 284.40816788543395d0*Lmu,& !1 ,x , x^2
      5432.878085716809d0 + 4617.753005620716d0*Lmu, &
	-2563.28806902233d0 - 2127.204749929736d0*Lmu, &
	-17.78991469587654d0 - 7.338839370492167d0*Lmu,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      -76.36755190995123d0 - 27.78047292742369d0*Lmu, &
	78.87763768206408d0 + 27.108672170009157d0*Lmu,&     !x/(1-x)Log^2(x), x Log^2(x)
      3443.1424448947796d0 + 2682.892435111155d0*Lmu, &
	-599.7983485750569d0 - 193.78350958785984d0*Lmu, &
	839.4963238597323d0 + 837.4027040063779d0*Lmu,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      544.0265746128981d0 + 488.0462345288209d0*Lmu, &
	1417.3624064262308d0 + 1497.1299554721447d0*Lmu, &
	-113.40180701234924d0 - 106.96326421760017d0*Lmu, &
	0.009338040302161874d0 + 0.00862378569953127d0*Lmu/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
    end if
  end subroutine Set_Coeff_q_qb
  
     !!!!!coefficient function q<-qp regular-part  
  subroutine Set_Coeff_q_qp(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qp=(/&
    0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
    0d0,0d0,&     !1/x, Log(x)/x
    0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
    0d0,0d0,0d0,& !1 ,x , x^2
    0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
    0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
    0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
    0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
      Coeff_q_qp=Coeff_q_qp+alpha*alpha*((/ &!!Exact
	0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
	688d0/81d0 - 208d0*Lmu/27d0 + 16d0*Lmu**2/9d0 - 32d0*zeta2/9d0, 0d0,&     !1/x, Log(x)/x
	8d0/3d0 - 8d0*Lmu/3d0 + 8d0*Lmu**2/3d0, -2d0/3d0 + 8d0*Lmu/3d0, 4d0/9d0,& !Log(x),Log^2(x),Log^3(x)
	32d0*Lmu/3d0 + 4d0*Lmu**2/3d0, -16d0*Lmu - 4d0*Lmu**2/3d0, 352d0*Lmu/27d0 - 16d0*Lmu**2/9d0,& !1 ,x , x^2
	0d0, -8d0*Lmu + 8d0*Lmu**2/3d0, -64d0*Lmu/9d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
	0d0, 8d0*Lmu/3d0,&     !x/(1-x)Log^2(x), x Log^2(x)
	0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
	0d0,0d0,0d0,0d0/)& !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
	+(/&!Approximate
	0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
	0d0,0d0,&     !1/x, Log(x)/x
	0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
	-603.924227035499d0, -4636.485211211568d0, -49.76555465398209d0,& !1 ,x , x^2
	-5287.52982020046d0, 2269.612280502602d0, -58.06494427244658d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
	119.75596348356056d0, -129.79997791500733d0,&     !x/(1-x)Log^2(x), x Log^2(x)
	-3369.724995744234d0, 427.8946185110229d0, -812.4665998422748d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
	-600.6972087253562d0, -1469.0061919285804d0, 92.73615445019001d0, -0.021528986881156446d0/)) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
    end if
  end subroutine Set_Coeff_q_qp
  
  !!! This function has been used during debuging
 subroutine CheckCoefficient(as,Nf,Lmu,z)
 real(dp)::Lmu,as,z,Nf,lz,l1z
 real(dp), dimension(1:23)::func
 real(dp), dimension(1:3)::func1
 
  lz=Log(z)
  l1z=Log(1d0-z)
  func=(/l1z,l1z**2,l1z**3,&
	    1d0/z,lz/z,&  
	    lz,lz**2,lz**3,&
	    1d0,z,z**2,& 
	    z/(1d0-z)*lz, z*lz, (z**2)*lz,& 
	    z/(1d0-z)*lz**2,z*lz**2,& 
	    (lz/(1d0-z)+1d0)*l1z,lz*l1z,z*lz*l1z,& 
	    (1d0-z)/z*l1z, (1d0-z)*l1z, ((1d0-z)**2)*l1z,(1d0-z)*l1z**2/)
     
  func1=(/1d0/(1d0-z),Log(1d0-z)/(1d0-z),Log(1d0-z)**2/(1d0-z)/)
 
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
  write(*,*) SUM(Coeff_g_g*func)+SUM(CoeffSing1_g_g*func1)
 end subroutine CheckCoefficient
