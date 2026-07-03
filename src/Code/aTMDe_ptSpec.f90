!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       arTeMiDe 3.05
!
!   Contains type definitions for kinematic specification of xSec points in cross-sections
!   and their updates.
!
!                               A.Vladimirov (29.06.2026)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module aTMDe_ptSpec
use aTMDe_numerics
implicit none
private

real(dp),parameter::zero=1.d-6 !!!! value of 0 to compare with


!!!!!! The specification of the Drell-Yan kinematic point.
type, public :: DYpoint
    private
    real(dp),public::s     !!!! DY Mandelshtam variable s
    real(dp),public::Q     !!!! photon's invariant mass
    real(dp),public::Q2    !!!! photon's invariant mass-squared (to avoid recomputing each time)
    real(dp),public::qT    !!!! photon's transverse momentum
    real(dp),public::y     !!!! photon's rapidity
    real(dp),public::expy  !!!! exponent of rapidity (to avoid recomputing each time)
    real(dp),public::tau   !!!! sqrt(q^+q^-/P^+P^-) = sqrt[(Q2+qT2)/s]

contains

end type

!!!!!! The specification of the SIDIS kinematic point.
type, public :: SIDISpoint
    private
    real(dp),public::sM2   !!!! s-M2 with SIDIS Mandelshtam variable s
    real(dp),public::Q     !!!! photon's invariant mass
    real(dp),public::Q2    !!!! photon's invariant mass-squared (to avoid recomputing each time)
    real(dp),public::pPerp !!!! produced hadron's transverse momentum
    real(dp),public::M2target!!!! mass-squared of the target hadron
    real(dp),public::M2product!!!! mass-squared of the product hadron
    real(dp),public::x    !!!! SIDIS variable x
    real(dp),public::y    !!!! SIDIS variable y
    real(dp),public::z    !!!! SIDIS variable z
    real(dp),public::gamma2    !!!! SIDIS variable gamma^2=(2xM/Q)^2

contains

end type

public::updateQ_DY,updateQT_DY,updateY_DY,coshY,yFromXF

public::updateQ_SIDIS,updatePperp_SIDIS,updateX_SIDIS,updateZ_SIDIS,eps_SIDIS,qT_SIDIS,ComputeX1Z1qT,pJacobian,QfromY,XfromY

interface DYpoint
    procedure :: DYpoint_constructor
end interface DYpoint

interface SIDISpoint
    procedure :: SIDISpoint_constructor
end interface SIDISpoint

contains
!!!!!!--------------------------------------------------------------------------------------------
!!!!!!-----------------------------------DYpoint related functions--------------------------------
!!!!!!--------------------------------------------------------------------------------------------

!!! Constructor for type DYpoint
!!! qT,Q,y are transverse momentum, invariant mass, and rapidity of the photon
!!! s is the Mandelshtam variable for DY
pure function DYpoint_constructor(qT,s,Q,y) result(this)
type(DYpoint)::this
real(dp), intent(in)::qT,s,Q,y

this%s=s
this%qT=qT
this%Q=Q
this%Q2=Q*Q
this%y=y
this%expy=exp(y)
this%tau=sqrt((this%Q2+qT*qT)/s)
end function DYpoint_constructor

!!!!!---- it is important to use these functions not as procedures of pointDY, but as stand alone subroutines
!!!!!---- in this case, one directly specifies type-variable, and saves time of looking through ID's of classes.
!!!!!---- although the overhead is small, it could be costly in integrations (very-many calls of updates)

!!!!! update the variable Q in the DYpoint
pure subroutine updateQ_DY(this,Q)
type(DYpoint),intent(inout)::this
real(dp),intent(in)::Q
this%Q=Q
this%Q2=Q*Q
this%tau=sqrt((this%Q2+this%qT*this%qT)/this%s)
end subroutine updateQ_DY

!!!!! update the variable qT in the DYpoint
pure subroutine updateQT_DY(this,qT)
type(DYpoint),intent(inout)::this
real(dp),intent(in)::qT
this%qT=qT
this%tau=sqrt((this%Q2+qT*qT)/this%s)
end subroutine updateQT_DY

!!!!! update the variable y in the DYpoint
pure subroutine updateY_DY(this,y)
type(DYpoint),intent(inout)::this
real(dp),intent(in)::y
this%y=y
this%expy=exp(y)
end subroutine updateY_DY

!!!!! returns cosh(y) computed from stored expy
pure function coshY(this)
type(DYpoint),intent(in)::this
real(dp)::coshY
coshY=(this%expy+1._dp/this%expy)/2
end function coshY

!!! Computes y-variable from x_F, for given kinematic point
pure function yFromXF(xF,this)
type(DYpoint),intent(in)::this
real(dp),intent(in):: xF
real(dp):: yFromXF
  yFromXF=asinh(xF/2d0/this%tau)
end function yFromXF

!!!!!!--------------------------------------------------------------------------------------------
!!!!!!----------------------------------SIDISpoint related functions------------------------------
!!!!!!--------------------------------------------------------------------------------------------

!!! Constructor for type SIDISpoint
!!! pPerp,Q,x,z are transverse momentum, invariant mass, x, and z for SIDIS
!!! s is the Mandelshtam variable for SIDIS
!!! Mtarget and Mproduct are corresponding masses of hadrons (in GeV)
pure function SIDISpoint_constructor(pPerp,s,Q,x,z,Mtarget,Mproduct) result(this)
type(SIDISpoint)::this
real(dp), intent(in)::pPerp,s,Q,x,z,Mtarget,Mproduct

this%pPerp=pPerp
this%Q=Q
this%Q2=Q*Q
this%x=x
this%z=z
this%M2target=Mtarget*Mtarget
this%M2product=Mproduct*Mproduct
this%sM2=s-this%M2target
this%y=this%Q2/(this%x*this%sM2)
this%gamma2=(2*x*Mtarget/Q)**2

end function SIDISpoint_constructor

!!!!!---- it is important to use these functions not as procedures of pointSIDIS, but as stand alone subroutines
!!!!!---- in this case, one directly specifies type-variable, and saves time of looking through ID's of classes.
!!!!!---- although the overhead is small, it could be costly in integrations (very-many calls of updates)

!!!!! update the variable Q in the SIDISpoint
pure subroutine updateQ_SIDIS(this,Q)
type(SIDISpoint),intent(inout)::this
real(dp),intent(in)::Q
this%Q=Q
this%Q2=Q*Q
this%y=this%Q2/(this%x*this%sM2)
this%gamma2=(2*this%x/Q)**2*this%M2target
end subroutine updateQ_SIDIS


!!!!! update the variable Pperp in the SIDISpoint
pure subroutine updatePperp_SIDIS(this,pPerp)
type(SIDISpoint),intent(inout)::this
real(dp),intent(in)::pPerp
this%pPerp=pPerp
end subroutine updatePperp_SIDIS

!!!!! update the variable x in the SIDISpoint
pure subroutine updateX_SIDIS(this,x)
type(SIDISpoint),intent(inout)::this
real(dp),intent(in)::x
this%x=x
this%y=this%Q2/(this%x*this%sM2)
this%gamma2=(2*this%x/this%Q)**2*this%M2target
end subroutine updateX_SIDIS

!!!!! update the variable z in the SIDISpoint
pure subroutine updateZ_SIDIS(this,z)
type(SIDISpoint),intent(inout)::this
real(dp),intent(in)::z
this%z=z
end subroutine updateZ_SIDIS

!!!!! returns SIDIS variable eps=(1-y-gamma^2 y^2/4)(1-y+y^2/2+gamma^2y^2/4)
pure function eps_SIDIS(this)
type(SIDISpoint),intent(in)::this
real(dp)::eps_SIDIS

real(dp)::y2,dummy1,dummy2
  y2 =this%y*this%y
  dummy1=1._dp-this%y
  dummy2=this%gamma2*y2*0.25_dp

  eps_SIDIS=(dummy1-dummy2)/(dummy1+y2*0.5_dp+dummy2)
end function eps_SIDIS

!!!!! returns SIDIS variable qT=pPper/z sqrt((1+gamma^2)(1-gamma2 gammah2))
pure function qT_SIDIS(this)
type(SIDISpoint),intent(in)::this
real(dp)::qT_SIDIS

!!!! in the case of massless hadron the formula is much simpler
if(this%M2target>zero) then
    if(this%M2product>zero) then
        qT_SIDIS=this%pPerp/this%z*sqrt((1+this%gamma2)/(1-this%gamma2*this%M2product/(this%z*this%z*this%Q2)))
    else
        qT_SIDIS=this%pPerp/this%z*sqrt(1+this%gamma2)
    end if
else
    qT_SIDIS=this%pPerp/this%z
end if
end function qT_SIDIS

!!!!! returns values of variables qT, x1 and z1 used as inputs to the structure functions
pure subroutine ComputeX1Z1qT(this,x1,z1,qT)
type(SIDISpoint),intent(in)::this
real(dp),intent(out)::qT,x1,z1

qT=qT_SIDIS(this)
x1=this%x*(1-qT*qT/this%Q2)
z1=this%z
end subroutine ComputeX1Z1qT

!!!!! returns SIDIS expression sqrt(1-gamma2 gamma2_h-gamma2 gamma_perp^2), which is a part of the Jacobian of transformation to p_perp
pure function pJacobian(this)
type(SIDISpoint),intent(in)::this
real(dp)::pJacobian

!!!! in the case of massless hadron the formula is much simpler
if(this%M2target>zero) then
    pJacobian=sqrt(1-this%gamma2*(this%M2product+this%pPerp*this%pPerp)/(this%z*this%z*this%Q2))
else
    pJacobian=1._dp
end if
end function pJacobian

!!!! Compute X based on given y and other parameters of the point [Q,sM2]
pure function XfromY(y,this)
type(SIDISpoint),intent(in)::this
real(dp),intent(in)::y
real(dp)::XfromY
XfromY=this%Q2/y/this%sM2
end function XfromY

!!!! Compute Q based on given x and other parameters of the point [x,sM2]
pure function QfromY(y,this)
type(SIDISpoint),intent(in)::this
real(dp),intent(in)::y
real(dp)::QfromY
QfromY=Sqrt(this%sM2*this%x*y)
end function QfromY

end module aTMDe_ptSpec
