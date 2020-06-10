!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDX_DY
implicit none

character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/'
character(*),parameter::repFILE=prefix//'Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep'
character(*),parameter::constFILE=prefix//'constants-files/DY+SIDIS_nnlo/const-DY+SIDIS_NNPDF31+DSS_nnlo'
real*8,allocatable::NParray(:)
integer::numR

integer,dimension(1:3)::proc
real*8::s,yMin,yMax,Qmin,Qmax,ptCut,ptCut2
real*8,allocatable::X0(:),Xplus(:),Xminus(:),XplusS(:),XminusS(:),qTmin(:),qTmax(:)
real*8::binSize
real::t1,t2,t3

! character(*),parameter::suffixFILE='.dat'
! character(*),parameter::suffixDescription='# Specification: Bins covered by theory'
character(*),parameter::suffixFILE='-full.dat'
character(*),parameter::suffixDescription='# Specification: all bins (including bins outside of theory range)'
real*8::time1,time2,time3

integer::i


call artemide_Initialize(constFILE)

numR=artemide_NumOfReplicasInFile(repFILE)

call artemide_GetReplicaFromFile(repFILE,0,NParray)
!call artemide_SetNPparameters_TMDR(NParray(1:2))
call artemide_SetNPparameters(NParray)

proc=(/1,1,5/) !!!! pp DY

s=13000d0**2

ymin=-2.4d0
ymax=2.4d0
ptCut=25d0
ptCut2=20d0

Qmin=20.d0
Qmax=40.d0

call TMDX_DY_setProcess(proc)
call TMDX_DY_XSetup(s,(Qmax+Qmin)/2d0,(ymax+ymin)/2d0)
call TMDX_DY_SetCuts(.true.,ptCut,ptCut2,yMin,yMax)

!!!!!!! original binning

call cpu_time(time1)
write(*,*) '---------------------------Q=50-76---------------------------------------'
call cpu_time(time2)

!!!!!!! Mass range 50-76 GeV: 16 bins:
!!!!!!! 0., 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 37., 52., 85., 160., 250., 1000.

Qmin=50.d0
Qmax=76.d0

allocate(qTmin(1:16))
allocate(qTmax(1:16))
qTmin=real((/0.0001, 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 37., 52., 85., 160., 250./),8)
qTmax=real((/2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 37., 52., 85., 160., 250., 1000./),8)

! allocate(qTmin(1:8))
! allocate(qTmax(1:8))
! qTmin=real((/0.0001, 2., 4., 6., 8., 10., 12., 14./),8)
! qTmax=real((/2., 4., 6., 8., 10., 12., 14., 18./),8)

call ComputeXsecWithStatErr(qTmin,qTmax,X0,Xminus,Xplus)
call ComputeXsecWithScaleErr(qTmin,qTmax,X0,XminusS,XplusS)


OPEN(UNIT=55, FILE="/home/vla18041/WorkingFiles/TMD/Fit_Notes/Predictions/aTMDe(Q=50-76)"//trim(suffixFILE)&
                    , ACTION="write", STATUS="replace")
    write(55,*) '# artemide results'
    write(55,*) trim(suffixDescription)
    write(55,*) '# s=',s,' GeV'
    write(55,*) '# Q=[',QMin,',',Qmax,'] GeV'
    write(55,*) '# Y=[',ymin,',',ymax,']'
    write(55,*) '# etaCut=[',ymin,',',ymax,']'
    write(55,*) '# pTCut=[',pTCut,',',ptCut2,']'
    write(55,*) '# Xsec= dsigma/dqT [pb/GeV]'
    write(55,*) "#qT(min)  qT(max)     Xsec   +dX(NP-param) -dX(NP-param) +dX(scale) -dX(scale)"
    do i=1,size(qTmin)
        binSize=qTmax(i)-qTmin(i)
        !write(*,'("{",F8.1,",",F8.1,",",F12.6,",",F12.6,",",F12.6,"},")') qTmin(i),qTmax(i),X0(i),Xminus(i),Xplus(i)
        write(*,'(F8.1,F8.1,F12.6,F12.6,F12.6,F12.6,F12.6)') qTmin(i),qTmax(i),&
            X0(i)/binSize, (Xminus(i)-X0(i))/binSize,(Xplus(i)-X0(i))/binSize,&
            (XminusS(i)-X0(i))/binSize,(XplusS(i)-X0(i))/binSize
        write(55,'(F8.1,F8.1,F12.6,F12.6,F12.6,F12.6,F12.6)') qTmin(i),qTmax(i),&
            X0(i)/binSize, (Xminus(i)-X0(i))/binSize,(Xplus(i)-X0(i))/binSize,&
            (XminusS(i)-X0(i))/binSize,(XplusS(i)-X0(i))/binSize
    end do
CLOSE(55,STATUS='KEEP')

call cpu_time(time3)
write(*,'("calc.time=",F10.2,"s. ")') time3-time2
write(*,*) ' '

write(*,*) '---------------------------Q=76-106---------------------------------------'
call cpu_time(time2)

!!!!!!! Mass range 76-106 GeV: 37 bins:
!!!!!!! 0.,   1.,   2.,   3.,   4.,   5.,   6.,    7.,   8.,   9.,
!!!!!!! 10.,  11.,  12.,  13.,  14.,  16.,  18.,   20.,  22.,  25.,
!!!!!!! 28.,  32.,  37.,  43.,  52.,  65.,  85.,  120., 160., 190.,
!!!!!!! 220., 250., 300., 350., 400., 450., 500., 1000.

Qmin=76.d0
Qmax=106.d0

deallocate(qTmin,qTmax)

allocate(qTmin(1:37))
allocate(qTmax(1:37))
qTmin=real((/0.0001, 1.,   2.,   3.,   4.,   5.,   6.,    7.,   8.,   9.,&
10.,  11.,  12.,  13.,  14.,  16.,  18.,   20.,  22.,  25.,&
28.,  32.,  37.,  43.,  52.,  65.,  85.,  120., 160., 190.,&
220., 250., 300., 350., 400., 450., 500./), 8)
qTmax=real((/1.,   2.,   3.,   4.,   5.,   6.,    7.,   8.,   9.,&
10.,  11.,  12.,  13.,  14.,  16.,  18.,   20.,  22.,  25.,&
28.,  32.,  37.,  43.,  52.,  65.,  85.,  120., 160., 190.,&
220., 250., 300., 350., 400., 450., 500., 1000./), 8)

! allocate(qTmin(1:19))
! allocate(qTmax(1:19))
! qTmin=real((/0.0001, 1.,   2.,   3.,   4.,   5.,   6.,    7.,   8.,   9.,&
! 10.,  11.,  12.,  13.,  14.,  16.,  18.,   20.,  22./), 8)
! qTmax=real((/1.,   2.,   3.,   4.,   5.,   6.,    7.,   8.,   9.,&
! 10.,  11.,  12.,  13.,  14.,  16.,  18.,   20.,  22.,  25./), 8)


call ComputeXsecWithStatErr(qTmin,qTmax,X0,Xminus,Xplus)
call ComputeXsecWithScaleErr(qTmin,qTmax,X0,XminusS,XplusS)
OPEN(UNIT=55, FILE="/home/vla18041/WorkingFiles/TMD/Fit_Notes/Predictions/aTMDe(Q=76-106)"//trim(suffixFILE)&
                    , ACTION="write", STATUS="replace")
    write(55,*) '# artemide results'
    write(55,*) trim(suffixDescription)
    write(55,*) '# s=',s,' GeV'
    write(55,*) '# Q=[',QMin,',',Qmax,'] GeV'
    write(55,*) '# Y=[',ymin,',',ymax,']'
    write(55,*) '# etaCut=[',ymin,',',ymax,']'
    write(55,*) '# pTCut=[',pTCut,',',ptCut2,']'
    write(55,*) '# Xsec= dsigma/dqT [pb/GeV]'
    write(55,*) "#qT(min)  qT(max)     Xsec   +dX(NP-param) -dX(NP-param) +dX(scale) -dX(scale)"
    do i=1,size(qTmin)
        binSize=qTmax(i)-qTmin(i)
        !write(*,'("{",F8.1,",",F8.1,",",F12.6,",",F12.6,",",F12.6,"},")') qTmin(i),qTmax(i),X0(i),Xminus(i),Xplus(i)
        write(*,'(F8.1,F8.1,F12.6,F12.6,F12.6,F12.6,F12.6)') qTmin(i),qTmax(i),&
            X0(i)/binSize, (Xminus(i)-X0(i))/binSize,(Xplus(i)-X0(i))/binSize,&
            (XminusS(i)-X0(i))/binSize,(XplusS(i)-X0(i))/binSize
        write(55,'(F8.1,F8.1,F12.6,F12.6,F12.6,F12.6,F12.6)') qTmin(i),qTmax(i),&
            X0(i)/binSize, (Xminus(i)-X0(i))/binSize,(Xplus(i)-X0(i))/binSize,&
            (XminusS(i)-X0(i))/binSize,(XplusS(i)-X0(i))/binSize
    end do
CLOSE(55,STATUS='KEEP')

call cpu_time(time3)
write(*,'("calc.time=",F10.2,"s. ")') time3-time2
write(*,*) ' '

write(*,*) '---------------------------Q=106-170---------------------------------------'
call cpu_time(time2)

!!!!!!! Mass range 106-170 GeV: 16 bins: 
!!!!!!! 0., 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 37., 52., 85., 160., 250., 1000.

Qmin=106.d0
Qmax=170.d0

deallocate(qTmin,qTmax)

allocate(qTmin(1:16))
allocate(qTmax(1:16))
qTmin=real((/0.0001, 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 37., 52., 85., 160., 250./),8)
qTmax=real((/2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 37., 52., 85., 160., 250., 1000./),8)

! allocate(qTmin(1:11))
! allocate(qTmax(1:11))
! qTmin=real((/0.0001, 2., 4., 6., 8., 10., 12., 14., 18., 22., 28./),8)
! qTmax=real((/2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 37./),8)

call ComputeXsecWithStatErr(qTmin,qTmax,X0,Xminus,Xplus)
call ComputeXsecWithScaleErr(qTmin,qTmax,X0,XminusS,XplusS)


OPEN(UNIT=55, FILE="/home/vla18041/WorkingFiles/TMD/Fit_Notes/Predictions/aTMDe(Q=106-170)"//trim(suffixFILE)&
                    , ACTION="write", STATUS="replace")
    write(55,*) '# artemide results'
    write(55,*) trim(suffixDescription)
    write(55,*) '# s=',s,' GeV'
    write(55,*) '# Q=[',QMin,',',Qmax,'] GeV'
    write(55,*) '# Y=[',ymin,',',ymax,']'
    write(55,*) '# etaCut=[',ymin,',',ymax,']'
    write(55,*) '# pTCut=[',pTCut,',',ptCut2,']'
    write(55,*) '# Xsec= dsigma/dqT [pb/GeV]'
    write(55,*) "#qT(min)  qT(max)     Xsec   +dX(NP-param) -dX(NP-param) +dX(scale) -dX(scale)"
    do i=1,size(qTmin)
        binSize=qTmax(i)-qTmin(i)
        !write(*,'("{",F8.1,",",F8.1,",",F12.6,",",F12.6,",",F12.6,"},")') qTmin(i),qTmax(i),X0(i),Xminus(i),Xplus(i)
        write(*,'(F8.1,F8.1,F12.6,F12.6,F12.6,F12.6,F12.6)') qTmin(i),qTmax(i),&
            X0(i)/binSize, (Xminus(i)-X0(i))/binSize,(Xplus(i)-X0(i))/binSize,&
            (XminusS(i)-X0(i))/binSize,(XplusS(i)-X0(i))/binSize
        write(55,'(F8.1,F8.1,F12.6,F12.6,F12.6,F12.6,F12.6)') qTmin(i),qTmax(i),&
            X0(i)/binSize, (Xminus(i)-X0(i))/binSize,(Xplus(i)-X0(i))/binSize,&
            (XminusS(i)-X0(i))/binSize,(XplusS(i)-X0(i))/binSize
    end do
CLOSE(55,STATUS='KEEP')

call cpu_time(time3)
write(*,'("calc.time=",F10.2,"s. ")') time3-time2
write(*,*) ' '

write(*,*) '---------------------------Q=170-350---------------------------------------'
call cpu_time(time2)

!!!!! Mass ranges 170-350 and 350-1000: 13 bins:
!!!!! 0., 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 52., 160., 1000.

Qmin=170.d0
Qmax=350.d0

deallocate(qTmin,qTmax)

allocate(qTmin(1:13))
allocate(qTmax(1:13))
qTmin=real((/0.0001, 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 52., 160./), 8)
qTmax=real((/2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 52., 160., 1000./), 8)

! allocate(qTmin(1:11))
! allocate(qTmax(1:11))
! qTmin=real((/0.0001, 2., 4., 6., 8., 10., 12., 14., 18., 22., 28./), 8)
! qTmax=real((/2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 52./), 8)

call ComputeXsecWithStatErr(qTmin,qTmax,X0,Xminus,Xplus)
call ComputeXsecWithScaleErr(qTmin,qTmax,X0,XminusS,XplusS)
OPEN(UNIT=55, FILE="/home/vla18041/WorkingFiles/TMD/Fit_Notes/Predictions/aTMDe(Q=170-350)"//trim(suffixFILE)&
                    , ACTION="write", STATUS="replace")
    write(55,*) '# artemide results'
    write(55,*) trim(suffixDescription)
    write(55,*) '# s=',s,' GeV'
    write(55,*) '# Q=[',QMin,',',Qmax,'] GeV'
    write(55,*) '# Y=[',ymin,',',ymax,']'
    write(55,*) '# etaCut=[',ymin,',',ymax,']'
    write(55,*) '# pTCut=[',pTCut,',',ptCut2,']'
    write(55,*) '# Xsec= dsigma/dqT [pb/GeV]'
    write(55,*) "#qT(min)  qT(max)     Xsec   +dX(NP-param) -dX(NP-param) +dX(scale) -dX(scale)"
    do i=1,size(qTmin)
        binSize=qTmax(i)-qTmin(i)
        !write(*,'("{",F8.1,",",F8.1,",",F12.6,",",F12.6,",",F12.6,"},")') qTmin(i),qTmax(i),X0(i),Xminus(i),Xplus(i)
        write(*,'(F8.1,F8.1,F12.6,F12.6,F12.6,F12.6,F12.6)') qTmin(i),qTmax(i),&
            X0(i)/binSize, (Xminus(i)-X0(i))/binSize,(Xplus(i)-X0(i))/binSize,&
            (XminusS(i)-X0(i))/binSize,(XplusS(i)-X0(i))/binSize
        write(55,'(F8.1,F8.1,F12.6,F12.6,F12.6,F12.6,F12.6)') qTmin(i),qTmax(i),&
            X0(i)/binSize, (Xminus(i)-X0(i))/binSize,(Xplus(i)-X0(i))/binSize,&
            (XminusS(i)-X0(i))/binSize,(XplusS(i)-X0(i))/binSize
    end do
CLOSE(55,STATUS='KEEP')


call cpu_time(time3)
write(*,'("calc.time=",F10.2,"s. ")') time3-time2
write(*,*) ' '

write(*,*) '---------------------------Q=350-1000---------------------------------------'
call cpu_time(time2)

!!!!! Mass ranges 170-350 and 350-1000: 13 bins:
!!!!! 0., 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 52., 160., 1000.

Qmin=350.d0
Qmax=1000.d0

deallocate(qTmin,qTmax)

allocate(qTmin(1:13))
allocate(qTmax(1:13))
qTmin=real((/0.0001, 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 52., 160./), 8)
qTmax=real((/2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 52., 160., 1000./), 8)

! allocate(qTmin(1:12))
! allocate(qTmax(1:12))
! qTmin=real((/0.0001, 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 52./), 8)
! qTmax=real((/2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 52., 160./), 8)

call ComputeXsecWithStatErr(qTmin,qTmax,X0,Xminus,Xplus)
call ComputeXsecWithScaleErr(qTmin,qTmax,X0,XminusS,XplusS)
OPEN(UNIT=55, FILE="/home/vla18041/WorkingFiles/TMD/Fit_Notes/Predictions/aTMDe(Q=350-1000)"//trim(suffixFILE)&
                    , ACTION="write", STATUS="replace")
    write(55,*) '# artemide results'
    write(55,*) trim(suffixDescription)
    write(55,*) '# s=',s,' GeV'
    write(55,*) '# Q=[',QMin,',',Qmax,'] GeV'
    write(55,*) '# Y=[',ymin,',',ymax,']'
    write(55,*) '# etaCut=[',ymin,',',ymax,']'
    write(55,*) '# pTCut=[',pTCut,',',ptCut2,']'
    write(55,*) '# Xsec= dsigma/dqT [pb/GeV]'
    write(55,*) "#qT(min)  qT(max)     Xsec   +dX(NP-param) -dX(NP-param) +dX(scale) -dX(scale)"
    do i=1,size(qTmin)
        binSize=qTmax(i)-qTmin(i)
        !write(*,'("{",F8.1,",",F8.1,",",F12.6,",",F12.6,",",F12.6,"},")') qTmin(i),qTmax(i),X0(i),Xminus(i),Xplus(i)
        write(*,'(F8.1,F8.1,F12.6,F12.6,F12.6,F12.6,F12.6)') qTmin(i),qTmax(i),&
            X0(i)/binSize, (Xminus(i)-X0(i))/binSize,(Xplus(i)-X0(i))/binSize,&
            (XminusS(i)-X0(i))/binSize,(XplusS(i)-X0(i))/binSize
        write(55,'(F8.1,F8.1,F12.6,F12.6,F12.6,F12.6,F12.6)') qTmin(i),qTmax(i),&
            X0(i)/binSize, (Xminus(i)-X0(i))/binSize,(Xplus(i)-X0(i))/binSize,&
            (XminusS(i)-X0(i))/binSize,(XplusS(i)-X0(i))/binSize
    end do
CLOSE(55,STATUS='KEEP')

call cpu_time(time3)
write(*,'("calc.time=",F10.2,"s. ")') time3-time2
write(*,*) ' '

write(*,'("TOTAL calc.time=",F10.2,"s. ")') time3-time1
write(*,*) ' '

contains 

!!! compute statistical error due to replica distribution
subroutine ComputeXsecWithStatErr(qTmin_in,qTmax_in,Xout,XminusOut,XplusOut)
  real*8,intent(in)::qTmin_in(:),qTmax_in(:)
  real*8,allocatable,intent(out)::XOut(:),XminusOut(:),XplusOut(:)
  real*8,allocatable::mean(:),deviation(:),X(:),Err(:)
  
  integer::l,i,j
  
  call artemide_SetScaleVariations(1d0,1d0,1d0,1d0)
  
  l=size(qTmin_in)
  allocate(XOut(1:l))
  allocate(XminusOut(1:l))
  allocate(XplusOut(1:l))
  allocate(mean(1:l))
  allocate(deviation(1:l))
  allocate(X(1:l))
  allocate(Err(1:l))
  
  do i=1,l
    mean(i)=0d0
    deviation(i)=0d0
  end do
  
  do j=1,numR
    call artemide_GetReplicaFromFile(repFILE,j,NParray)
    call artemide_SetNPparameters(NParray)
      
    call CalcXsec_DY_PTint_Qint_Yint(X,qTmin_in,qTmax_in,Qmin,Qmax,yMin,yMax)
    
    mean=mean+X
    deviation=deviation+X**2
    
  end do

  do i=1,l
    Xout(i)=mean(i)/numR
    Err(i)=Sqrt(deviation(i)/numR - Xout(i)**2)
    XminusOut(i)=Xout(i)-Err(i)
    XplusOut(i)=Xout(i)+Err(i)
  end do
  
end subroutine ComputeXsecWithStatErr

!!! compute scle-variation error
subroutine ComputeXsecWithScaleErr(qTmin_in,qTmax_in,XOut,XminusOut,XplusOut)
  real*8,intent(in)::qTmin_in(:),qTmax_in(:)
  real*8,allocatable,intent(out)::XOut(:),XminusOut(:),XplusOut(:)
  real*8,allocatable::X(:),X1(:),X2(:),X3(:),X4(:)
  real*8::r
  integer::l,i,j
  
  call artemide_GetReplicaFromFile(repFILE,0,NParray) 
  call artemide_SetNPparameters(NParray)

  
  l=size(qTmin_in)
  allocate(XOut(1:l))
  allocate(XminusOut(1:l))
  allocate(XplusOut(1:l))
  
  allocate(X(1:l))
  allocate(X1(1:l))
  allocate(X2(1:l))
  allocate(X3(1:l))
  allocate(X4(1:l))
  
  call artemide_SetScaleVariations(1d0,1d0,1d0,1d0)
  call CalcXsec_DY_PTint_Qint_Yint(X,qTmin_in,qTmax_in,Qmin,Qmax,yMin,yMax)

  r=2d0
  
  call artemide_SetScaleVariations(1d0,1d0/r,1d0,1d0)
  call CalcXsec_DY_PTint_Qint_Yint(X1,qTmin_in,qTmax_in,Qmin,Qmax,yMin,yMax)
  call artemide_SetScaleVariations(1d0,1d0*r,1d0,1d0)
  call CalcXsec_DY_PTint_Qint_Yint(X2,qTmin_in,qTmax_in,Qmin,Qmax,yMin,yMax)
  call artemide_SetScaleVariations(1d0,1d0,1d0,1d0/r)
  call CalcXsec_DY_PTint_Qint_Yint(X3,qTmin_in,qTmax_in,Qmin,Qmax,yMin,yMax)
  call artemide_SetScaleVariations(1d0,1d0,1d0,1d0*r)
  call CalcXsec_DY_PTint_Qint_Yint(X4,qTmin_in,qTmax_in,Qmin,Qmax,yMin,yMax)
  
  do i=1,l
    Xout(i)=X(i)
    XminusOut(i)=min(X(i),X1(i),X2(i),X3(i),X4(i))
    XplusOut(i)=max(X(i),X1(i),X2(i),X3(i),X4(i))
  end do

  
end subroutine ComputeXsecWithScaleErr

end program example
