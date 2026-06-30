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

public::updateQ_DY,updateQT_DY,updateY_DY,coshY,yFromXF

interface DYpoint
    procedure :: DYpoint_constructor
end interface DYpoint

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

end module aTMDe_ptSpec



