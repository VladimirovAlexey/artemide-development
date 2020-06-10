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
real*8::s,yMin,yMax,ptCut,ptCut2
real*8,allocatable::X0(:),Xplus(:),Xminus(:),XplusS(:),XminusS(:)
real*8,allocatable::Q0(:),qT0(:)
real*8,allocatable::Qlist(:,:),qTlist(:,:),yList(:,:),CutPList(:,:),sList(:)
integer,allocatable::procList(:,:),NumList(:)
logical,allocatable::inCutList(:)
real*8::binSize,delta
real::t1,t2,t3

character(*),parameter::suffixFILE='-fine.dat'
character(*),parameter::suffixDescription='# Specification: even finer pt-bins and finer Q-bins'
real*8::time1,time2,time3

integer::i,j,k, iQ,iqT,modNum


proc=(/1,1,5/) !!!! pp DY

s=13000d0**2

ymin=-2.4d0
ymax=2.4d0
ptCut=25d0
ptCut2=20d0

!------------------------------------
!----- create lists of kinematics
iQ=10
allocate(Q0(1:iQ))
Q0=(/50d0,76d0,81d0,86d0,91d0,96d0,106d0,170d0,350d0,1000d0/)

! iQ=19
! allocate(Q0(1:iQ))
! Q0=(/50d0,63d0,76d0,79d0,81d0,84d0,86d0,89d0,91d0,93d0,96d0,101d0,106d0,138d0,170d0,260d0,350d0,675d0,1000d0/)

iqT=38
allocate(qT0(1:iqT))
qT0=(/0.0001d0, 1d0,   2d0,   3d0,   4d0,   5d0,   6d0,    7d0,   8d0,   9d0,&
 10d0,  11d0,  12d0,  13d0,  14d0,  16d0,  18d0,   20d0,  22d0,  25d0,&
 28d0,  32d0,  37d0,  43d0,  52d0,  65d0,  85d0,  120d0, 160d0, 190d0,&
 220d0, 250d0, 300d0, 350d0, 400d0, 450d0, 500d0,1000d0/)
! iqT=91
! allocate(qT0(1:iqT))
! qT0=(/0.0001d0, 0.5d0, 1d0,  1.5d0, 2d0, 2.5d0,  3d0, 3.5d0,  4d0, 4.5d0, 5d0, &
!  5.5d0, 6d0, 6.5d0, 7d0, 7.5d0, 8d0, 8.5d0, 9d0,9.5d0,&
!  10d0,  11d0,  12d0,  13d0,  14d0,  15d0,  16d0,   17d0,  18d0,  19d0,&
!  20d0, 21d0, 22d0, 23d0, 24d0, 25d0, 26d0, 27d0, 28d0, 29d0, 30d0, &
!  32d0, 34d0, 36d0, 38d0, 40d0, 42d0, 44d0,46d0,48d0,50d0,52d0, 54d0, 56d0, 58d0, 60d0, 62d0, 64d0,66d0,68d0,&
!  70d0,75d0,80d0,85d0,90d0,95d0,100d0,110d0,120d0,130d0,140d0,150d0,160d0,180d0,200d0,&
!  220d0,240d0,260d0,280d0,300d0,350d0,400d0,450d0,500d0,550d0,600d0,650d0,700d0,800d0,900d0,1000d0/)

 ! count number of bins
 k=0
 delta=0.25d0
 do i=1,iQ-1
 do j=1,iqT-1
 if((qT0(j+1)+qT0(j))<delta*(Q0(i+1)+Q0(i))) k=k+1
 end do
 end do
 
 write(*,*) '-----------------------------Total number of bins:', k,'---------------------'
 allocate(slist(1:k))
 allocate(Qlist(1:k,1:2))
 allocate(qTList(1:k,1:2))
 allocate(yList(1:k,1:2))
 allocate(CutPList(1:k,1:4))
 allocate(procList(1:k,1:3))
 allocate(NumList(1:k))
 allocate(inCutList(1:k))
 !!! create lists
 k=1
 do i=1,iQ-1
 do j=1,iqT-1
 if((qT0(j+1)+qT0(j))<delta*(Q0(i+1)+Q0(i))) then
  procList(k,1:3)=proc
  sList(k)=s
  Qlist(k,1:2)=(/Q0(i),Q0(i+1)/)
  qTlist(k,1:2)=(/qT0(j),qT0(j+1)/)
  ylist(k,1:2)=(/yMin,yMax/)
  inCutList(k)=.true.
  CutPList(k,1:4)=(/ptCut,ptCut2,yMin,yMax/)
  modNum=Int((qT0(j+1)-qT0(j))/5)
  if(modNum>4) then
    if(mod(modNum,2)==1) then
     modNum=modNum+1
    end if
    NumList(k)=modNum
  else
    NumList(k)=4
  end if  
  k=k+1
 end if
 end do
 end do
 

call artemide_Initialize(constFILE)

numR=artemide_NumOfReplicasInFile(repFILE)

call artemide_GetReplicaFromFile(repFILE,0,NParray)
!call artemide_SetNPparameters_TMDR(NParray(1:2))
call artemide_SetNPparameters(NParray)


call TMDX_DY_setProcess(proc)
call TMDX_DY_XSetup(s,91d0,(ymax+ymin)/2d0)
call TMDX_DY_SetCuts(.true.,ptCut,ptCut2,yMin,yMax)

!!!!!!! original binning

call cpu_time(time1)
write(*,*) '---------------------------COMPUTATION START---------------------------------------'
call cpu_time(time2)

call ComputeXsecWithStatErr(X0,Xminus,Xplus)
call cpu_time(time3)
write(*,'("(statistical err) calc.time=",F10.2,"s. ")') time3-time2
write(*,*) ' '
call cpu_time(time2)
call ComputeXsecWithScaleErr(X0,XminusS,XplusS)
call cpu_time(time3)
write(*,'("(scale err) calc.time=",F10.2,"s. ")') time3-time2


OPEN(UNIT=55, FILE="/home/vla18041/WorkingFiles/TMD/Fit_Notes/Predictions/aTMDe"//trim(suffixFILE)&
                    , ACTION="write", STATUS="replace")
    write(55,*) '# artemide results'
    write(55,*) trim(suffixDescription)
    write(55,*) '# s=',s,' GeV'
    write(55,*) '# Y=[',ymin,',',ymax,']'
    write(55,*) '# etaCut=[',ymin,',',ymax,']'
    write(55,*) '# pTCut=[',pTCut,',',ptCut2,']'
    write(55,*) '# Xsec= dsigma/dqT [pb/GeV]'
    write(55,*) "#Q(min) Q(max) qT(min)  qT(max)     Xsec   +dX(NP-param) -dX(NP-param) +dX(scale) -dX(scale)"
    do i=1,size(sList)
        binSize=qTlist(i,2)-qTlist(i,1)
        write(*,'(F8.1,F8.1,F8.1,F8.1,F12.6,F12.6,F12.6,F12.6,F12.6)') Qlist(i,1),Qlist(i,2), qTList(i,1),qTList(i,2),&
            X0(i)/binSize, (Xminus(i)-X0(i))/binSize,(Xplus(i)-X0(i))/binSize,&
            (XminusS(i)-X0(i))/binSize,(XplusS(i)-X0(i))/binSize
        write(55,'(F8.1,F8.1,F8.1,F8.1,F12.6,F12.6,F12.6,F12.6,F12.6)') Qlist(i,1),Qlist(i,2), qTList(i,1),qTList(i,2),&
            X0(i)/binSize, (Xminus(i)-X0(i))/binSize,(Xplus(i)-X0(i))/binSize,&
            (XminusS(i)-X0(i))/binSize,(XplusS(i)-X0(i))/binSize
    end do
CLOSE(55,STATUS='KEEP')

call cpu_time(time3)

write(*,'("TOTAL calc.time=",F10.2,"s. ")') time3-time1
write(*,*) ' '

contains 

!!! compute statistical error due to replica distribution
subroutine ComputeXsecWithStatErr(Xout,XminusOut,XplusOut)
  real*8,allocatable,intent(out)::XOut(:),XminusOut(:),XplusOut(:)
  real*8,allocatable::mean(:),deviation(:),X(:),Err(:)
  
  integer::l,i,j
  
  call artemide_SetScaleVariations(1d0,1d0,1d0,1d0)
  
  l=size(sList)
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
      
    call xSec_DY_List(X,procList,sList,qTList,Qlist,yList,inCutList,CutPList,NumList)
    
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
subroutine ComputeXsecWithScaleErr(XOut,XminusOut,XplusOut)
  real*8,allocatable,intent(out)::XOut(:),XminusOut(:),XplusOut(:)
  real*8,allocatable::X(:),X1(:),X2(:),X3(:),X4(:)
  real*8::r
  integer::l,i,j
  
  call artemide_GetReplicaFromFile(repFILE,0,NParray) 
  call artemide_SetNPparameters(NParray)

  
  l=size(sList)
  allocate(XOut(1:l))
  allocate(XminusOut(1:l))
  allocate(XplusOut(1:l))
  
  allocate(X(1:l))
  allocate(X1(1:l))
  allocate(X2(1:l))
  allocate(X3(1:l))
  allocate(X4(1:l))
  
  call artemide_SetScaleVariations(1d0,1d0,1d0,1d0)
  call xSec_DY_List(X,procList,sList,qTList,Qlist,yList,inCutList,CutPList,NumList)

  r=2d0
  
  call artemide_SetScaleVariations(1d0,1d0/r,1d0,1d0)
  call xSec_DY_List(X1,procList,sList,qTList,Qlist,yList,inCutList,CutPList,NumList)
  call artemide_SetScaleVariations(1d0,1d0*r,1d0,1d0)
  call xSec_DY_List(X2,procList,sList,qTList,Qlist,yList,inCutList,CutPList,NumList)
  call artemide_SetScaleVariations(1d0,1d0,1d0,1d0/r)
  call xSec_DY_List(X3,procList,sList,qTList,Qlist,yList,inCutList,CutPList,NumList)
  call artemide_SetScaleVariations(1d0,1d0,1d0,1d0*r)
  call xSec_DY_List(X4,procList,sList,qTList,Qlist,yList,inCutList,CutPList,NumList)
  
  do i=1,l
    Xout(i)=X(i)
    XminusOut(i)=min(X(i),X1(i),X2(i),X3(i),X4(i))
    XplusOut(i)=max(X(i),X1(i),X2(i),X3(i),X4(i))
  end do

  
end subroutine ComputeXsecWithScaleErr

end program example
