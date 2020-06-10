program test_cuts
use LeptonCutsDY
implicit none

integer::i,j,k
real*8::Q,qT,y

  write(*,*) '--------------'
  do i=0, 10
  write(*,*) 0.0d0+i*0.2,CutFactor4(0.4d0,91d0,0.0d0+i*0.2,(/20d0,20d0,-2.4d0,2.4d0/)),&
      CutFactorA(0.4d0,91d0,0.0d0+i*0.2,(/20d0,20d0,-2.4d0,2.4d0/))
  end do
  
    write(*,*) '--------------'
  do i=0, 10
  write(*,*) 2.0d0+i*0.2,CutFactor4(0.4d0,91d0,2.0d0+i*0.2,(/20d0,20d0,2.d0,4.5d0/)),&
      CutFactorA(0.4d0,91d0,2.0d0+i*0.2,(/20d0,20d0,2d0,4.5d0/))
  end do
  
  do i=0,10
  do j=0,10
  do k=0,10
   qT=1d0*i
   Q=91d0-(j-5)*5d0
   y=0.0d0+k*0.2
   write(*,*) qT,Q,y,CutFactorA(qT,Q,y,(/20d0,20d0,-2.4d0,2.4d0/))/CutFactor4(qT,Q,y,(/20d0,20d0,-2.4d0,2.4d0/))
  end do
  end do
  end do
  
end program