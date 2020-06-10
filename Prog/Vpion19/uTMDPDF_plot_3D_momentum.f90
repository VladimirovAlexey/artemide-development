program example
use aTMDe_control
use uTMDPDF
use TMDs_inKT
implicit none

character(*),parameter::constFILE=&
'/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Prog/Vpion19/const-piDY+kT'

real*8,allocatable::NParray(:)
integer::numR,numB,i,j
real*8,allocatable::central(:),mean(:),deviation(:),kT(:)
real*8::kMax,step,mu,dd,x,TMD(-5:5)
integer::f,k
!real*8,parameter::xValues(1:7)=(/.99d0,10d0**(-0.5d0),10d0**(-1d0),10d0**(-1.5d0),10d0**(-2d0),10d0**(-2.5d0),10d0**(-3d0)/)
real*8,parameter::xValues(1:13)=(/&
1.d0, 0.923077d0, 0.846154d0, 0.769231d0, 0.692308d0, 0.615385d0, 0.538462d0,&
0.461538d0, 0.384615d0, 0.307692d0, 0.230769d0, 0.153846d0, 0.0769231d0/)

call artemide_Initialize(constFILE)
numR=100

call uTMDPDF_SetLambdaNP(0,.true.,.false.)


numB=40
kMax=3d0
mu=4d0

allocate(kT(1:numB))
do i=1,20
  kT(i)=0.005d0+i/20d0
end do
do i=21,numB
  kT(i)=kMax+(kMax-1d0)*(i-numB)/(numB-20)
end do

allocate(central(1:numB))
allocate(mean(1:numB))
allocate(deviation(1:numB))

write(*,*) ' '
write(*,*) '----------------------------D-QUARK----------------------------------------------'
write(*,*) ' '

f=1

do k=1,size(xValues)
x=xValues(k)

  do i=1,numB
    TMD=uTMDPDF_kT_5(x,kT(i),2)
    central(i)=TMD(f)
    mean(i)=0d0
    deviation(i)=0d0
  end do

  do j=1,numR
    
    call uTMDPDF_SetLambdaNP(j,.true.,.false.)
    do i=1,numB
      TMD=uTMDPDF_kT_5(x,kT(i),2)
      mean(i)=mean(i)+TMD(f)
      deviation(i)=deviation(i)+TMD(f)**2
    end do  
    
  end do

  do i=1,numB
    mean(i)=mean(i)/numR
    deviation(i)=Sqrt(deviation(i)/numR-mean(i)**2)
  end do


  do i=1,numB
    write(*,"('{'F10.7,',',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") x,kT(i),x*mean(i),&
	  x*(mean(i)-deviation(i)),x*(mean(i)+deviation(i))
  end do

end do

write(*,*) ' '
write(*,*) '----------------------------U-QUARK----------------------------------------------'
write(*,*) ' '

f=2

do k=1,size(xValues)
x=xValues(k)

  do i=1,numB
    TMD=uTMDPDF_kT_5(x,kT(i),2)
    central(i)=TMD(f)
    mean(i)=0d0
    deviation(i)=0d0
  end do

  do j=1,numR
    !call artemide_GetReplicaFromFile(repFILE,j,NParray)
    !call uTMDPDF_SetLambdaNP(NParray(3:9),.true.,.false.)
    call uTMDPDF_SetLambdaNP(j,.true.,.false.)
    do i=1,numB
      TMD=uTMDPDF_kT_5(x,kT(i),2)
      mean(i)=mean(i)+TMD(f)
      deviation(i)=deviation(i)+TMD(f)**2
    end do  
    
  end do

  do i=1,numB
    mean(i)=mean(i)/numR
    deviation(i)=Sqrt(deviation(i)/numR-mean(i)**2)
  end do


  do i=1,numB
    write(*,"('{'F10.7,',',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") x,kT(i),x*mean(i),&
	  x*(mean(i)-deviation(i)),x*(mean(i)+deviation(i))
  end do

end do

write(*,*) ' '
write(*,*) '----------------------------dBAR-QUARK----------------------------------------------'
write(*,*) ' '

f=-1

do k=1,size(xValues)
x=xValues(k)

  do i=1,numB
    TMD=uTMDPDF_kT_5(x,kT(i),2)
    central(i)=TMD(f)
    mean(i)=0d0
    deviation(i)=0d0
  end do

  do j=1,numR
    !call artemide_GetReplicaFromFile(repFILE,j,NParray)
    !call uTMDPDF_SetLambdaNP(NParray(3:9),.true.,.false.)
    call uTMDPDF_SetLambdaNP(j,.true.,.false.)
    do i=1,numB
      TMD=uTMDPDF_kT_5(x,kT(i),2)
      mean(i)=mean(i)+TMD(f)
      deviation(i)=deviation(i)+TMD(f)**2
    end do  
    
  end do

  do i=1,numB
    mean(i)=mean(i)/numR
    deviation(i)=Sqrt(deviation(i)/numR-mean(i)**2)
  end do


  do i=1,numB
    write(*,"('{'F10.7,',',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") x,kT(i),x*mean(i),&
	  x*(mean(i)-deviation(i)),x*(mean(i)+deviation(i))
  end do

end do

write(*,*) ' '
write(*,*) '----------------------------UBAR-QUARK----------------------------------------------'
write(*,*) ' '

f=-2

do k=1,size(xValues)
x=xValues(k)

  do i=1,numB
    TMD=uTMDPDF_kT_5(x,kT(i),2)
    central(i)=TMD(f)
    mean(i)=0d0
    deviation(i)=0d0
  end do

  do j=1,numR
    !call artemide_GetReplicaFromFile(repFILE,j,NParray)
    !call uTMDPDF_SetLambdaNP(NParray(3:9),.true.,.false.)
    call uTMDPDF_SetLambdaNP(j,.true.,.false.)
    do i=1,numB
      TMD=uTMDPDF_kT_5(x,kT(i),2)
      mean(i)=mean(i)+TMD(f)
      deviation(i)=deviation(i)+TMD(f)**2
    end do  
    
  end do

  do i=1,numB
    mean(i)=mean(i)/numR
    deviation(i)=Sqrt(deviation(i)/numR-mean(i)**2)
  end do


  do i=1,numB
    write(*,"('{'F10.7,',',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") x,kT(i),x*mean(i),&
	  x*(mean(i)-deviation(i)),x*(mean(i)+deviation(i))
  end do

end do

end program example
