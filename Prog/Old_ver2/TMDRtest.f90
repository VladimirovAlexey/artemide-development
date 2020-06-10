program example
use TMDR
implicit none

real*8,dimension(1:2):: A,B
real*8::bT,Q1,Q2,Q3 ,x1,x2,x3,dd
integer::i
 A=(/3.98d0,0.2d0/)
 B=(/3.5217d0,0.0293d0/)
  
  call TMDR_Initialize('const-test')
 
 bT=0.24d0
 Q1=4.125d0
 Q2=25d0
 Q3=100d0
 
 A=(/6d0,0.0d0/)
 call TMDR_setNPparameters(A)
  
!   write(*,*) DNP(3d0,100d0,1),DNP(40d0,100d0,1),DNP(400d0,100d0,1)
!  
!  A=(/5d0,0.0d0/)
!   call TMDR_setNPparameters(A)
!   
!   write(*,*) DNP(3d0,100d0,1),DNP(40d0,100d0,1),DNP(400d0,100d0,1)
!   
!    A=(/4d0,0.0d0/)
!   call TMDR_setNPparameters(A)
!   
!   write(*,*) DNP(3d0,100d0,1),DNP(40d0,100d0,1),DNP(400d0,100d0,1)
!   
!    A=(/3d0,0.0d0/)
!   call TMDR_setNPparameters(A)
!   
!   write(*,*) DNP(3d0,100d0,1),DNP(40d0,100d0,1),DNP(400d0,100d0,1)
!   
!    A=(/2d0,0.0d0/)
!   call TMDR_setNPparameters(A)
!   
!   write(*,*) DNP(3d0,100d0,1),DNP(40d0,100d0,1),DNP(400d0,100d0,1)
!   
!    A=(/1d0,-0.1d0/)
!   call TMDR_setNPparameters(A)
!   
!   write(*,*) DNP(3d0,100d0,1),DNP(40d0,100d0,1),DNP(400d0,100d0,1)
!   
!   stop
  
!   x1=TMDR_Rzeta(bT,Q1,Q1**2,1)
!   x2=TMDR_Rzeta(bT,Q2,Q2**2,1)
!   x3=TMDR_Rzeta(bT,Q3,Q3**2,1)
!   write(*,*) x1,x2,x3
!   
!   bT=0.5d0
! 
!   x1=TMDR_Rzeta(bT,Q1,Q1**2,1)
!   x2=TMDR_Rzeta(bT,Q2,Q2**2,1)
!   x3=TMDR_Rzeta(bT,Q3,Q3**2,1)
!   write(*,*) x1,x2,x3
!   
!   bT=1.0d0
! 
!   x1=TMDR_Rzeta(bT,Q1,Q1**2,1)
!   x2=TMDR_Rzeta(bT,Q2,Q2**2,1)
!   x3=TMDR_Rzeta(bT,Q3,Q3**2,1)
!   write(*,*) x1,x2,x3
!   
!   bT=2.0d0
! 
!   x1=TMDR_Rzeta(bT,Q1,Q1**2,1)
!   x2=TMDR_Rzeta(bT,Q2,Q2**2,1)
!   x3=TMDR_Rzeta(bT,Q3,Q3**2,1)
!   write(*,*) x1,x2,x3
!   write(*,*) '------------------------------------------------'
! !   
!   dd=0.1d0
!   bt=1d0
!   x1=zetaMUpert(Q2,bt,1)
!   x2=zetaMuResum(Q2,bt,1)
!   x3=zetaMUresum3(Q2,bt,1)
!   write(*,*) 'bt=',bt,' : ',x1,x2,x3
!   
!   bt=0.3d0
!   x1=zetaMUpert(Q2,bt,1)
!   x2=zetaMuResum(Q2,bt,1)
!   x3=zetaMUresum3(Q2,bt,1)
!   write(*,*) 'bt=',bt,' : ',x1,x2,x3
!   
!   bt=0.6d0
!   x1=zetaMUpert(Q2,bt,1)
!   x2=zetaMuResum(Q2,bt,1)
!   x3=zetaMUresum(Q2,bt,1)
!   write(*,*) 'bt=',bt,' : ',x1,x2,x3
!   
!   bt=1.0d0
!   x1=zetaMUpert(Q2,bt,1)
!   x2=zetaMuResum(Q2,bt,1)
!   x3=zetaMUresum3(Q2,bt,1)
!   write(*,*) 'bt=',bt,' : ',x1,x2,x3
!   
!   bt=1.5d0
!   x1=zetaMUpert(Q2,bt,1)
!   x2=zetaMuResum(Q2,bt,1)
!   x3=zetaMUresum3(Q2,bt,1)
!   write(*,*) 'bt=',bt,' : ',x1,x2,x3

!--------------------------------------------
  call TMDR_setNPparameters(B)
  
  Q1=103d0
  bt=25.d0
!   x1= zetaSL(Q1,bt,1)
!   write(*,*) x1
  
  Q1=120d0
  do i=0,10
    bt=0.1+i*5d0
    x1=TMDR_R_toSL(bt,Q1,Q1**2,1)
    x2=TMDR_R_toSL(bt,Q1,zetaNP(Q1,bt,1),1)
    x3=TMDR_Rzeta(bt,Q1,Q1**2,91d0,1)
    !write(*,*)bt,Q1,x1,x2,x1/x2,log(zetaNP(Q1,bt,1)/zetaSL(Q1,bt,1)),TMDR_Rzeta(bt,Q1,Q1**2,91d0,1)
    write(*,*)bt,Q1,x1,x2,x1/x2,x3
  end do

end program example