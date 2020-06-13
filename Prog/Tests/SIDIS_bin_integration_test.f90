!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program tests the bin integration for SIDIS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDX_SIDIS
implicit none

real*8::time1,time2

real*8,allocatable::xSec(:),s(:),pT(:,:),z(:,:),x(:,:),Q(:,:),Cuts(:,:),masses(:,:)
logical,allocatable::doCut(:)
integer,allocatable::proc_list(:,:)

integer::i,iMax


call artemide_Initialize('const-SIDIS_test',prefix='Prog/Tests/const-files/')
call artemide_SetNPparameters_uTMDPDF((/0.253434d0, 9.04351d0, 346.999d0, 2.47992d0, -5.69988d0, 0.d0, 0.d0/))

call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.253434d0, 9.04351d0, 346.999d0, 2.47992d0, -5.69988d0, 0.d0, 0.d0/))
call artemide_SetNPparameters_uTMDFF((/0.264d0,0.479d0,0.459d0,0.539d0/))

iMax=5

allocate(xSec(1:iMax),s(1:iMax),pT(1:iMax,1:2),z(1:iMax,1:2),x(1:iMax,1:2),Q(1:iMax,1:2),Cuts(1:iMax,1:4),masses(1:iMax,1:2))
allocate(doCut(1:iMax))
allocate(proc_list(1:iMax,1:3))


write(*,*) '----------------------------------------------------------------------------------------------'
write(*,*) '-------------------------------------Simple test ---------------------------------------------'
write(*,*) '---------------------- s=300, 4<Q<6, 0.3<z<0.5, 0.1<x<0.2 ------------------------------------'
write(*,*)



do i=1,iMax
    proc_list(i,1:3)=(/1,1,2001/) !!!! p->pi+ SIDIS
    s(i)=300d0
    Q(i,1:2)=(/4d0,6d0/)
    z(i,1:2)=(/0.3d0,0.5d0/)
    x(i,1:2)=(/0.1d0,0.2d0/)
    pt(i,1:2)=(/0.1d0*i,0.1d0*(i+1)/)
    doCut(i)=.false.
    cuts(i,1:4)=(/0.1d0,0.9d0,25d0,10000d0/)
    masses(i,1:2)=(/0.938d0,0.12d0/)
end do

call cpu_time(time1)
call xSec_SIDIS_List(xsec,proc_list,s,pT,z,x,Q,doCut,Cuts,masses)
call cpu_time(time2)
write(*,'("|       pt         |   xSec       |")')
do i=1,iMax
  write(*,'("| [",F6.2,",",F6.2,"]  |",F12.8,"  |")') pt(i,1),pt(i,2),Xsec(i)/1000d0
end do
write(*,'("   TOTAL TIME:     ",F10.4)') time2-time1

write(*,*) '----------------------------------------------------------------------------------------------'
write(*,*) '-------------------------------------Large x-range -------------------------------------------'
write(*,*) '---------------------- s=300, 4<Q<6, 0.3<z<0.5, 0.01<x<1 -------------------------------------'
write(*,*)

do i=1,iMax
    proc_list(i,1:3)=(/1,1,2001/) !!!! p->pi+ SIDIS
    s(i)=300d0
    Q(i,1:2)=(/4d0,6d0/)
    z(i,1:2)=(/0.3d0,0.5d0/)
    x(i,1:2)=(/0.01d0,1d0/)
    pt(i,1:2)=(/0.1d0*i,0.1d0*(i+1)/)
    doCut(i)=.false.
    cuts(i,1:4)=(/0.1d0,0.9d0,25d0,10000d0/)
    masses(i,1:2)=(/0.938d0,0.12d0/)
end do

call cpu_time(time1)
call xSec_SIDIS_List(xsec,proc_list,s,pT,z,x,Q,doCut,Cuts,masses)
call cpu_time(time2)

write(*,'("|       pt         |   xSec       |")')
do i=1,iMax
  write(*,'("| [",F6.2,",",F6.2,"]  |",F12.8,"  |")') pt(i,1),pt(i,2),Xsec(i)/1000d0
end do
write(*,'("   TOTAL TIME:     ",F10.4)') time2-time1

write(*,*) '----------------------------------------------------------------------------------------------'
write(*,*) '-------------------------------------Large z-range -------------------------------------------'
write(*,*) '---------------------- s=300, 4<Q<6, 0.2<z<0.9, 0.1<x<0.2 ------------------------------------'
write(*,*)

do i=1,iMax
    proc_list(i,1:3)=(/1,1,2001/) !!!! p->pi+ SIDIS
    s(i)=300d0
    Q(i,1:2)=(/4d0,6d0/)
    z(i,1:2)=(/0.2d0,0.9d0/)
    x(i,1:2)=(/0.1d0,0.2d0/)
    pt(i,1:2)=(/0.1d0*i,0.1d0*(i+1)/)
    doCut(i)=.false.
    cuts(i,1:4)=(/0.1d0,0.9d0,25d0,10000d0/)
    masses(i,1:2)=(/0.938d0,0.12d0/)
end do

call cpu_time(time1)
call xSec_SIDIS_List(xsec,proc_list,s,pT,z,x,Q,doCut,Cuts,masses)
call cpu_time(time2)

write(*,'("|       pt         |   xSec       |")')
do i=1,iMax
  write(*,'("| [",F6.2,",",F6.2,"]  |",F12.8,"  |")') pt(i,1),pt(i,2),Xsec(i)/1000d0
end do
write(*,'("   TOTAL TIME:     ",F10.4)') time2-time1

write(*,*) '----------------------------------------------------------------------------------------------'
write(*,*) '-------------------------------------Large Q-range -------------------------------------------'
write(*,*) '---------------------- s=300, 3<Q<25, 0.3<z<0.5, 0.1<x<0.2 -----------------------------------'
write(*,*)

do i=1,iMax
    proc_list(i,1:3)=(/1,1,2001/) !!!! p->pi+ SIDIS
    s(i)=300d0
    Q(i,1:2)=(/3d0,25d0/)
    z(i,1:2)=(/0.3d0,0.5d0/)
    x(i,1:2)=(/0.1d0,0.2d0/)
    pt(i,1:2)=(/0.1d0*i,0.1d0*(i+1)/)
    doCut(i)=.false.
    cuts(i,1:4)=(/0.1d0,0.9d0,25d0,10000d0/)
    masses(i,1:2)=(/0.938d0,0.12d0/)
end do

call cpu_time(time1)
call xSec_SIDIS_List(xsec,proc_list,s,pT,z,x,Q,doCut,Cuts,masses)
call cpu_time(time2)

write(*,'("|       pt         |   xSec       |")')
do i=1,iMax
  write(*,'("| [",F6.2,",",F6.2,"]  |",F12.8,"  |")') pt(i,1),pt(i,2),Xsec(i)/1000d0
end do
write(*,'("   TOTAL TIME:     ",F10.4)') time2-time1

end program example
