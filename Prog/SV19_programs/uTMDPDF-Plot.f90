!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Program that make a plot of unpolarized TMD F1
!!  for given x, mu, zeta in the range b in (0,..)
!!  saves into the file as (b, F(-6:6_)
!!
!!			AV: 16.09.2018
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use uTMDPDF
use aTMDe_control
implicit none
 character(*),parameter::name="Fx01_2_2.dat"
 integer :: i,j,f,k
 real*8 :: time1, time2
 integer,parameter::numB=40
 integer,parameter::numR=100
 real*8,allocatable::NParray(:)
 
 real*8,parameter :: fixB(0:numB)=(/&
  0.0039d0, 0.0046d0, 0.0055d0, 0.0066d0, 0.0078d0, 0.0093d0, 0.011d0, 0.0131d0, &
  0.0156d0, 0.0186d0, 0.0221d0, 0.0263d0, 0.0312d0, 0.0372d0, 0.0442d0, 0.0526d0, &
  0.0625d0, 0.0743d0, 0.0884d0, 0.1051d0, 0.125d0, 0.1487d0, 0.1768d0, 0.2102d0, 0.25d0, &
  0.2973d0, 0.3536d0, 0.4204d0, 0.5d0, 0.5946d0, 0.7071d0, 0.8409d0, 1.d0, 1.1892d0, &
  1.4142d0, 1.6818d0, 2.d0, 2.3784d0, 2.8284d0, 3.3636d0, 4.d0/)
 
 real*8 :: mu, zeta,x,b,TMD(-5:5)
 real*8 :: replicas(1:numR,0:numB,1:6),center(0:numB,1:6),mean(0:numB,1:6),deviation(0:numB,1:6)
 real*8 :: a1,a2,a3,a4,a5,a6,a7,a8,a9
 
 call uTMDPDF_Initialize('const-DY_TEST',&
      '/misc/data2/braun/vla18041/arTeMiDe_Repository/Constants-files/')
 
 
 call cpu_time(time1)
 
 
 call artemide_GetReplicaFromFile("/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19-replicas/DY_HERA20_nnlo.rep"&
	,0,NParray)
	
  call uTMDPDF_SetLambdaNP(NParray(3:9),.false.,.false.)
  x=0.5d0
  f=-1
 
 !!! central
 do i=0,numB
  mean(i,1:6)=(/0d0,0d0,0d0,0d0,0d0,0d0/)
  deviation(i,1:6)=(/0d0,0d0,0d0,0d0,0d0,0d0/)
  b=fixB(i)!real(i)/8d0+0.001d0
  TMD=uTMDPDF_lowScale5(x,b,1)
  center(i,1:6)=(/TMD(-3),TMD(-2),TMD(-1),TMD(1),TMD(2),TMD(3)/)
  !write(*,"('{',F10.7,',',F10.7,'},')") b,TMD(1)
 end do
 
!  ! replicas set by hand
!  do j=1,numR
!   !call uTMDPDF_SetPDFReplica(j)
! 
!   do i=0,numB        
!     b=fixB(i)!real(i)/8d0+0.001d0
!     TMD=uTMDPDF_lowScale5(x,b,1)
!     replicas(j,i)=TMD(f)
!     mean(i)=mean(i)+TMD(f)
!     deviation(i)=deviation(i)+TMD(f)**2
!     !write(*,*) b,TMD(1),mean(i),deviation(i)
!   end do
!  end do

!OPEN(UNIT=51, FILE="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/FIGURES_DY+SIDIS_2019/PDF-replicas/NNPDF_dist+.txt",&
!  OPEN(UNIT=51, FILE="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/FIGURES_DY+SIDIS_2019/PDF-replicas/NNPDF_stat.txt",& 
!       ACTION="read", STATUS="old")
 !! replicas
 do j=1,numR
!   read(51,*) k,a1,a2,a3,a4,a5,a6,a7,a8,a9  
  !call uTMDPDF_SetPDFReplica(k)  
  !call uTMDPDF_SetLambdaNP((/a3,a4,a5,a6,a7,a8,a9/),.false.,.false.)
!   call uTMDPDF_SetLambdaNP(j,.false.,.false.)
    call artemide_GetReplicaFromFile(&
    "/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19-replicas/DY_HERA20_nnlo.rep"&
	,j,NParray)
  call uTMDPDF_SetLambdaNP(NParray(3:9),.false.,.false.)
  do i=0,numB        
    b=fixB(i)!real(i)/8d0+0.001d0
    TMD=uTMDPDF_lowScale5(x,b,1)
    replicas(j,i,1:6)=(/TMD(-3),TMD(-2),TMD(-1),TMD(1),TMD(2),TMD(3)/)
    mean(i,1:6)=mean(i,1:6)+(/TMD(-3),TMD(-2),TMD(-1),TMD(1),TMD(2),TMD(3)/)
    deviation(i,1:6)=deviation(i,1:6)+(/TMD(-3)**2,TMD(-2)**2,TMD(-1)**2,TMD(1)**2,TMD(2)**2,TMD(3)**2/)
    !write(*,*) b,TMD(1),mean(i),deviation(i)
  end do
 end do
 
!  CLOSE (51, STATUS='KEEP')
 
 do i=0,numB
  mean(i,1:6)=mean(i,1:6)/numR
  deviation(i,1:6)=Sqrt(deviation(i,1:6)/numR-mean(i,1:6)**2)
 end do
 write(*,*) '{{'
 !write(*,*) '------------------------------anti s-----------------------------------------------'
 do f=1,6
 do i=0,numB
  b=fixB(i)!real(i)/8d0+0.001d0    
  if(i==numB) then
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'}}')") b,mean(i,f),mean(i,f)-deviation(i,f),mean(i,f)+deviation(i,f)
  else
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") b,mean(i,f),mean(i,f)-deviation(i,f),mean(i,f)+deviation(i,f)
  end if
 end do
 if(f==6) then
  write(*,*) '}'
 else
  write(*,*) ',{'
 end if
 end do
 
 call cpu_time(time2)

 write(*,*) 'timing=',time2-time1
 
end program example