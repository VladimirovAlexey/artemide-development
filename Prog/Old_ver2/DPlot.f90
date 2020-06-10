!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Program that make a plot rapidity anomalous dimension
!!  for given mu, in the range b in (0,..)
!!
!!	to make it working make DNP of TMDR public
!!
!!			AV: 16.09.2018
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use TMDR
implicit none
 character(*),parameter::name4="D_4.dat"
 character(*),parameter::name91="D_91.dat"
 integer :: j,k
 integer,parameter::num=100
 integer,parameter::numRep=95
 real*8 :: time1, time2
 real*8 :: mu,b
 real*8,dimension(1:num)::d0,dMean,deltaD
 real*8,dimension(1:num,1:numRep)::drep
 
 
 
 call cpu_time(time1)
 
 call TMDR_Initialize('const-HERAPDF20_NNLO','/home/vla18041/LinkData2/arTeMiDe_Repository/')
 
 mu=4d0
 
 !! mean value
 call TMDR_SetNPParameters(0)
 write(*,*) DNP(mu,0.1d0,1)
 
 do j=1,num
  b=real(j)/20d0
  d0(j)=DNP(mu,b,1)
 end do
 
 ! values for repilcas
 do k=1,numRep
 call TMDR_SetNPParameters(k)
 do j=1,num
  b=real(j)/20d0
  drep(j,k)=DNP(mu,b,1)
 end do
 end do
 
 !!! evaluate the mean values
 do j=1,num
  dMean(j)=sum(drep(j,1:numRep))/numRep
 end do
 write(*,*) "CHECK THE NUMBERS", (d0-dMean)/Abs(d0)
 
  !!! evaluate the variance
 do j=1,num
  deltaD(j)=sqrt(sum((drep(j,1:numRep)-d0(j))**2)/(numRep-1))
 end do
 
 open(7, file=name4, status="replace", action="write")
 do j=1,num
  b=real(j)/20d0
  write(7,*) b, d0(j),d0(j)-deltaD(j),d0(j)+deltaD(j)
 end do
close(7, STATUS='KEEP')
 
 
 call cpu_time(time2)

 write(*,*) 'timing=',time2-time1
 
 Write(*,*) '----------------------------------------------------------------------------------'
 
 mu=91d0
 
 !! mean value
 call TMDR_SetNPParameters(0)
 write(*,*) DNP(mu,0.1d0,1)
 
 do j=1,num
  b=real(j)/20d0
  d0(j)=DNP(mu,b,1)
 end do
 
 ! values for repilcas
 do k=1,numRep
 call TMDR_SetNPParameters(k)
 do j=1,num
  b=real(j)/20d0
  drep(j,k)=DNP(mu,b,1)
 end do
 end do
 
 !!! evaluate the mean values
 do j=1,num
  dMean(j)=sum(drep(j,1:numRep))/numRep
 end do
 write(*,*) "CHECK THE NUMBERS", (d0-dMean)/Abs(d0)
 
  !!! evaluate the variance
 do j=1,num
  deltaD(j)=sqrt(sum((drep(j,1:numRep)-d0(j))**2)/(numRep-1))
 end do
 
 open(7, file=name91, status="replace", action="write")
 do j=1,num
  b=real(j)/20d0
  write(7,*) b, d0(j),d0(j)-deltaD(j),d0(j)+deltaD(j)
 end do
close(7, STATUS='KEEP')
 
 
 call cpu_time(time2)

 write(*,*) 'timing=',time2-time1
 
end program example