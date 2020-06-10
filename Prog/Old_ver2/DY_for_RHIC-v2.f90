!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Calculation of DY cross-section for CMS.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program xSec_DY
use aTMDe_control
use TMDX_DY
implicit none

real*8 :: time1, time2
integer :: dummyINT,j
real*8 :: dummyREAL
logical::exist

integer:: process,ss
real*8:: Qmin,Qmax,s,ymin,ymax,ptCut,ptCut2
integer::s1,s2,s3,s4 !!! TMD sizes
real*8,allocatable, dimension(:) :: pt, xSec1
real*8,allocatable, dimension(:) :: xSec11!!!xSec for variation of c2+
real*8,allocatable, dimension(:) :: xSec12!!!xSec for variation of c2-
real*8,allocatable, dimension(:) :: xSec13!!!xSec for variation of c4+
real*8,allocatable, dimension(:) :: xSec14!!!xSec for variation of c4-


real*8,dimension(1:12)::ptBINS=(/0.1d0, 1.25d0, 2.5d0, 3.75d0, 5d0, 7.5d0, 10d0, 12.5d0, 15d0, 17.5d0, 20d0,25d0/)

call cpu_time(time1)

process=1!!

!This is from the artemide ver.1.4
!   call TMDX_DY_Initialize("NNLO")
!   call TMDX_DY_SetNPParameters((/3.3235d0,0.0380d0,0.2204d0, 7.0808d0,351.7950d0, 2.4632d0,-3.8334d0,0.0001d0, 0.0000d0/))

!This is for artemide ver.2.00
  call artemide_Initialize('constants-file','Models/BSV19.bFIT/')
!   call artemide_Initialize('const-DYfit18_LO','/home/vla18041/LinkData2/arTeMiDe_Repository/')
  call artemide_SetReplica_uTMDPDF(0)
  call artemide_SetReplica_TMDR(0)

     
  s=510d0**2

  ymin=-1.1d0
  ymax=1.1d0
  ptCut=25d0
  ptCut2=25d0
     
   call TMDX_DY_SetProcess(process)
   !call SetCuts(.true.,ptCut,ptCut2,yMin,yMax)
   call TMDX_DY_SetCuts(.true.,ptCut,ptCut2,yMin,yMax)
   call TMDX_DY_XSetup(s,91d0,0d0)
   
   !---------------------------------------------------------CENTRAL VALUE
   
  Qmin=73d0
  Qmax=126d0
  
   do j=1,11
   if(0.2d0*Qmax<ptBINS(j)) exit
   end do
   s1=j-1
   
   allocate(pt(1:j))
   allocate(xSec1(1:j-1),xSec11(1:j-1),xSec12(1:j-1),xSec13(1:j-1),xSec14(1:j-1))
   pt=ptBINS(1:j)
   ss=j-1
   
!    do j=1,ss
!     call CalcXsec_DY_PTint_Qint_Yint(xSec1(j),pt(j),pt(j+1),Qmin,Qmax,ymin,ymax)
!     write(*,*) pt(j:j+1),'-------',xSec1(j)
!    end do
   call CalcXsec_DY_PTint_Qint_Yint(xSec1,pt,Qmin,Qmax,ymin,ymax)
   xSec1=xSec1/(pt(2:ss+1)-pt(1:ss))
   deallocate(pt)
!-----------------------------------------------------------------
   
!---------------------------------------------------------Variation of c2+
  
!   call TMDX_DY_SetScaleVariations(1d0,2d0,1d0,1d0)
  call artemide_SetScaleVariations(1d0,2d0,1d0,1d0)
   
  Qmin=73d0
  Qmax=126d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec11,ptBINS(1:s1+1),Qmin,Qmax,ymin,ymax)
   xSec11=xSec11/(ptBins(2:s1+1)-ptBins(1:s1))
   
   !-----------------------------------------------------------------
   !---------------------------------------------------------Variation of c2-
  
!   call TMDX_DY_SetScaleVariations(1d0,0.5d0,1d0,1d0)
  call artemide_SetScaleVariations(1d0,0.5d0,1d0,1d0)
   
  Qmin=73d0
  Qmax=126d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec12,ptBINS(1:s1+1),Qmin,Qmax,ymin,ymax)
   xSec12=xSec12/(ptBins(2:s1+1)-ptBins(1:s1))
   
  !-----------------------------------------------------------------
  !---------------------------------------------------------Variation of c4+
  
!   call TMDX_DY_SetScaleVariations(1d0,1d0,1d0,2d0)
  call artemide_SetScaleVariations(1d0,1d0,1d0,2d0)
   
  Qmin=73d0
  Qmax=126d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec13,ptBINS(1:s1+1),Qmin,Qmax,ymin,ymax)
   xSec13=xSec13/(ptBins(2:s1+1)-ptBins(1:s1))
   
   !-----------------------------------------------------------------
   !---------------------------------------------------------Variation of c4-
  
!   call TMDX_DY_SetScaleVariations(1d0,1d0,1d0,0.5d0)
  call artemide_SetScaleVariations(1d0,1d0,1d0,0.5d0)
   
  Qmin=73d0
  Qmax=126d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec14,ptBINS(1:s1+1),Qmin,Qmax,ymin,ymax)
   xSec14=xSec14/(ptBins(2:s1+1)-ptBins(1:s1))
   
   
   !-----------------------------------------------------------------
   
   write(*,*) ' sqrt(s) = ',Sqrt(s)
   write(*,*) ' y in (',ymin,',',ymax,')'
   write(*,*) ' pt(cut1) = ',ptCut
   write(*,*) ' pt(cut2) = ',ptCut2
   
   write(*,*) ' '
   
   
   write(*,*) '73<Q< 126:  number of TMD bins ',size(xSec1)
   write(*,*) 'pt(min),   pt(max),    xSec,    err-,     err+'
   do j=1,size(xSec1)
   write(*,'(F6.2,"     ", F6.2,"     ", ES12.3E3,"     ", ES12.3E3,"     ", ES12.3E3)') ptBINS(j),ptBINS(j+1),xSec1(j),&
    min(xSec1(j),xSec11(j),xSec12(j),xSec13(j),xSec14(j)),&
    max(xSec1(j),xSec11(j),xSec12(j),xSec13(j),xSec14(j))
   end do
   
   write(*,*) "  "
   
call cpu_time(time2)
! 
write(*,*) 'Calculation time=',time2-time1
! 
   
 end program xSec_DY