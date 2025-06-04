!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Simple test-code for Snowlake. Checks compilation integrity and simple evaluation        !!
!!                                          A.Vladimirov 01.04.2024                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program SnowflakeTEST
use Snowflake
use Snowflake_Model

implicit none

real*8::mu0,mu1,Q(0:20),dd(0:20),xx(0:20)
integer::i,ff(0:20)

!call  SnowFlake_Initialize("Snowflake.ini")
call  SnowFlake_Initialize("prog/TEST.ini")

call SnowFlake_Model_Initialize()

call SetNPparameters(&
    (/0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0/))

mu0=1.d0
mu1=25.d0
call ComputeEvolution(mu0,mu1,alpha,U1=Tu,D1=Td,S1=Ts,U2=dTu,D2=dTd,S2=dTs,G1=Tp,G2=Tm,inputQ="T",inputG="T")
!call ComputeEvolution(mu0,mu1,alpha,U1=Tu,U2=dTu,inputQ="T",inputG="T")
!
! write(*,*) Tu(0.2d0,0.1d0),Tu(0.2d0,-0.1d0),Tu(0.3d0,0.1d0)
! write(*,*) Td(0.2d0,0.1d0),Td(0.2d0,-0.1d0),Td(0.3d0,0.1d0)
! write(*,*) Ts(0.2d0,0.1d0),Ts(0.2d0,-0.1d0),Ts(0.3d0,0.1d0)
! write(*,*) dTu(0.2d0,0.1d0),dTu(0.2d0,-0.1d0),dTu(0.3d0,0.1d0)
! write(*,*) dTd(0.2d0,0.1d0),dTd(0.2d0,-0.1d0),dTd(0.3d0,0.1d0)
! write(*,*) dTs(0.2d0,0.1d0),dTs(0.2d0,-0.1d0),dTs(0.3d0,0.1d0)
! write(*,*) dTpp(0.2d0,0.1d0),dTpp(0.2d0,-0.1d0),dTpp(0.3d0,0.1d0)
! write(*,*) dTmm(0.2d0,0.1d0),dTmm(0.2d0,-0.1d0),dTmm(0.3d0,0.1d0)



do i=0,20
xx(i)=0.3d0+i*0.01
Q(i)=2.d0+i
ff(i)=2
end do

write(*,*) '------'

call D2_list(dd,Q,ff)
do i=0,20
write(*,'("{",F5.2,",",F12.6,",",F12.6,"},")') Q(i),D2(Q(i),ff(i)),dd(i)
end do

write(*,*) '------'

call G2_list(dd,xx,Q,ff)
do i=0,20
write(*,'("{",F5.2,",",F12.6,",",F12.6,"},")') Q(i),G2(xx(i),Q(i),ff(i)),dd(i)
end do

end program SnowflakeTEST
