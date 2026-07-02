!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Part of code that contains hard coefficient functions for SIDIS-like x-Sections
!!                    is a part of artemide.TMDX_SIDIS_point
!!
!!            02.07.2026    created A.Vladimirov
!!
!!                                A.Vladimirov
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! hard coefficient taken from 1004.3653 up to 2-loop
!!! it takes global values of Q,order
!!! NOTE it uses Nf=3(fixed)
function HardCoefficientSIDIS(mu)
real(dp)::HardCoefficientSIDIS,mu,alpha,LQ!=Log[Q^2/mu^2]=-2Log[c1]

HardCoefficientSIDIS=1.d0
if(orderH_global>=1) then
    LQ=-2d0*LOG(c2_global)
    alpha=As(mu*c2_global)
    HardCoefficientSIDIS=HardCoefficientSIDIS+alpha*&
    (-16.946842488404727d0 + 8d0*LQ - 2.6666666666666665d0*LQ**2)
if(orderH_global>=2) then
    HardCoefficientSIDIS=HardCoefficientSIDIS+alpha**2*&
    (-116.50054911601637d0 + 46.190372772820254d0*LQ + 16.843858371984233d0*LQ**2&
    -13.333333333333334d0*LQ**3 + 3.5555555555555554d0*LQ**4)
if(orderH_global>=3) then
    HardCoefficientSIDIS=HardCoefficientSIDIS+alpha**3*&
    (-4820.715927678687 + 2492.274201933993*LQ + 44.19495641116441*LQ**2 &
    - 237.22228827339313*LQ**3 + 43.33848430014775*LQ**4 + 7.111111111111111*LQ**5 &
    -3.1604938271604937*LQ**6)
if(orderH_global>=4) then
    HardCoefficientSIDIS=HardCoefficientSIDIS+alpha**4*&
    (26391.759725461765 - 31391.21814540276*LQ + 16136.794429475773*LQ**2 &
    - 3922.584672164565*LQ**3 - 98.8284343286739*LQ**4 + 250.360398809412*LQ**5 &
    - 47.07700456533447*LQ**6 - 3.950617283950617*LQ**7 + 1.8436213991769548*LQ**8)
end if
end if
end if
end if
end function HardCoefficientSIDIS
