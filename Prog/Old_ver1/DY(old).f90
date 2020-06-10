!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Program that make DY cross-section (VERY OD)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program DY
use TMDX
implicit none

real*8 :: time1, time2
character(256)::line
integer :: dummyINT,j
real*8 :: dummyREAL

integer:: process, order,varyMOD,length,Xtype,qqq
real*8:: Q1,Q2,s, y1,y2
real*8:: pt,eta1,eta2,err1,err2
logical:: includeCuts,exist
real*8,allocatable, dimension(:) ::qtList1,qtList2,xSec,xSec1,xSec2,xSec3,xSec4,xSec5,xSec6,xSec7,xSec8

call cpu_time(time1)
  write(*,*) '---------------DY arTeMiDe-----------------------'
  Xtype=0
  write(*,"(A)",advance="no") '..read input'
    !!!!!!! --------------Read different parameters-------------------------
    OPEN(UNIT=51, FILE='bin/input', ACTION="read", STATUS="old")
    !!!!search for process
    do
    read(51,'(A)') line    
    if(line(1:3)=='*1 ') exit
    end do
    
    read(51,*) dummyINT
    if(dummyINT==1) then
     process=5
    else if(dummyINT==2) then
     process=7
    else
    write(*,*) 'ERROR IN INPUT LINE *1'
    stop
   end if
   
   !!!!search for s-entry
   do
    read(51,'(A)') line    
    if(line(1:3)=='*2 ') exit
   end do
    
    read(51,*) dummyREAL
    s=dummyREAL**2
   
   !!!!search for Q-entry
   do
    read(51,'(A)') line    
    if(line(1:3)=='*3 ') exit
   end do
    
    read(51,*) Q1,Q2
    if(Q2>0d0) then
     Xtype=Xtype+1   
    
    if(Q2<Q1) then
    dummyREAL=Q1
    Q1=Q2
    Q2=dummyREAL
    end if
    end if
    
   !!!!search for y-entry
   do
    read(51,'(A)') line    
    if(line(1:3)=='*4 ') exit
   end do
    
    read(51,*) y1,y2
    if(y2==0d0) then
     Xtype=Xtype+0
    else    
     Xtype=Xtype+10
     if(y1>y2) then
      dummyREAL=y1
      y1=y2
      y2=dummyREAL
     end if
    end if
    
    !!!!search for cuts
    do
     read(51,'(A)') line    
     if(line(1:3)=='*5 ') exit
    end do
    
    read(51,*) dummyINT
    if(dummyINT==1) then
     includeCuts=.true.
    else
     includeCuts=.false.
    end if
    
    if(includeCuts) then
     do
      read(51,'(A)') line    
      if(line(1:3)=='*5a ') exit
     end do
     read(51,*) pt
     do
      read(51,'(A)') line    
      if(line(1:3)=='*5b ') exit
     end do
     read(51,*) eta1,eta2
    end if
    
    !!!!search for order
   do
    read(51,'(A)') line    
    if(line(1:3)=='*6 ') exit
   end do
    
    read(51,*) order
    
        !!!!search for variation mode
   do
    read(51,'(A)') line    
    if(line(1:3)=='*7 ') exit
   end do
    
    read(51,*) varyMOD
    
   !!!!search for qt integration
   do
    read(51,'(A)') line    
    if(line(1:3)=='*8 ') exit
   end do
    read(51,*) qqq
     
     if(qqq>0) then
     Xtype=Xtype+100
     end if
     
   do
    read(51,'(A)') line    
    if(line(1:3)=='*8a') exit
   end do
    read(51,*) length
    
    allocate(qtList1(1:length))
    allocate(qtList2(1:length))
    allocate(xSec(1:length))
    allocate(xSec1(1:length))
    allocate(xSec2(1:length))
    allocate(xSec3(1:length))
    allocate(xSec4(1:length))
    allocate(xSec5(1:length))
    allocate(xSec6(1:length))
    allocate(xSec7(1:length))
    allocate(xSec8(1:length))
    
   do
    read(51,'(A)') line    
    if(line(1:3)=='*qT') exit
   end do
    
    do j=1,length
     if(qqq==0) then
      read(51,*) qtList1(j)
     else 
      read(51,*) qtList1(j),qtList2(j)
     end if
    end do
    
    CLOSE (51, STATUS='KEEP') 
    
    write(*,*) '   ..done.'
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SELECT CASE(order)
      CASE (1)
	call TMDX_Initialize("LO")
      CASE (2)
	call TMDX_Initialize("NLO")
      CASE (3)
	call TMDX_Initialize("NNLO")
      CASE DEFAULT
	write(*,*) 'UNKNOWN ORDER SWITCH TO NLO'
	call TMDX_Initialize("NLO")
     END SELECT
     
      write(*,*) 'sqrt(s)=',sqrt(s)
    if(Q2==0) then 
    write(*,*) 'Q=',Q1
    else
    write(*,*) 'Q=',Q1,' to ',Q2
    end if
    if(y2==0) then 
    write(*,*) 'y=',y1
    else
    write(*,*) 'y=',y1,' to ',y2
    end if
        SELECT CASE(order)
      CASE (1)
	write(*,*) 'Using LO'
      CASE (2)
	write(*,*) 'Using NLO'
      CASE (3)
	write(*,*) 'Using NNLO'
      CASE DEFAULT
	write(*,*) 'Using NLO'
     END SELECT
     
     if(includeCuts) then 
      call SetCuts(.true.,pt,eta1,eta2)
      write(*,*) 'Lepton cuts included with pt>',pt, eta1,'<eta<',eta2
     else
      call SetCuts(.false.,0d0,-1d0,1d0)
     end if
    
   call TMDX_XSetup(s,Q1,y1,process)
   call InitializeScaleVariations(1d0,1d0,1d0,1d0)
   
   
   !!!!THESE are parameters for MODEL1
     SELECT CASE(order)
      CASE (1)
	call TMDX_DY_SetNPParameters((/1.67d0,0.327d0,0.112d0,0.828d0/))
      CASE (2)
	call TMDX_DY_SetNPParameters((/2.43d0,0.106d0,0.183d0,0.328d0/))
      CASE (3)
	call TMDX_DY_SetNPParameters((/2.5208d0,0.0827d0,0.2446d0,0.3107d0/))
      CASE DEFAULT
	call TMDX_DY_SetNPParameters((/2.43d0,0.106d0,0.183d0,0.328d0/))
     END SELECT
   SELECT CASE(Xtype)
    CASE(0)
     write(*,*) '.. Evaluating dS/(dQ^2 dy dqT^2)'
     call CalculateXsec(xSec,qtList1)
    CASE(1)
     write(*,*) '.. Evaluating \INT dQ^2  dS/(dQ^2 dy dqT^2))'
     call CalculateXsec_Qint(xSec,qtList1,Q1,Q2)
    CASE(10)
     write(*,*) '.. Evaluating \INT dy  dS/(dQ^2 dy dqT^2))'
     call CalculateXsec_Yint(xSec,qtList1,y1,y2)
    CASE(100)
     write(*,*) '.. Evaluating \INT dqT^2  dS/(dQ^2 dy dqT^2))'
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
     stop
    CASE(11)
     write(*,*) '.. Evaluating \INT dQ^2 dy dS/(dQ^2 dy dqT^2))'
     call CalculateXsec_Qint_Yint(xSec,qtList1,y1,y2,Q1,Q2)
    CASE(101)
     write(*,*) '.. Evaluating \INT dQ^2 dqT^2 dS/(dQ^2 dy dqT^2))'
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
     stop
    CASE(110)
     write(*,*) '.. Evaluating 1/Delta qT \INT dqT^2 dy dS/(dQ^2 dy dqT^2))      (avarage over qT)'
     call CalculateXsec_PTint_Yint(xSec,qtList1,qtList2,y1,y2)
     xSec=xSec/(qtList2-qtList1)
    CASE(111)
     write(*,*) '.. Evaluating 1/Delta qT \INT dqT^2 dy dQ^2 dS/(dQ^2 dy dqT^2))      (avarage over qT)'
     call CalculateXsec_PTint_Qint_Yint(xSec,qtList1,qtList2,Q1,Q2,y1,y2)
     xSec=xSec/(qtList2-qtList1)
    CASE DEFAULT
     write(*,*) 'Unexpected combination of integrations....'
     stop
    END SELECT
    
    if(varyMOD>0) write(*,*) 'Evaluation of cross-section is done. Evaluation of error-band started.'
    
    If(varyMOD==1 .or. varyMOD==5 ) then
     call InitializeScaleVariations(0.5d0,1d0,1d0,1d0)
     SELECT CASE(Xtype)
    CASE(0)
     call CalculateXsec(xSec1,qtList1)
    CASE(1)
     call CalculateXsec_Qint(xSec1,qtList1,Q1,Q2)
    CASE(10)
     call CalculateXsec_Yint(xSec1,qtList1,y1,y2)
    CASE(100)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(11)
     call CalculateXsec_Qint_Yint(xSec1,qtList1,y1,y2,Q1,Q2)
    CASE(101)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(110)
     call CalculateXsec_PTint_Yint(xSec1,qtList1,qtList2,y1,y2)
     xSec1=xSec1/(qtList2-qtList1)
    CASE(111)
     call CalculateXsec_PTint_Qint_Yint(xSec1,qtList1,qtList2,Q1,Q2,y1,y2)
     xSec1=xSec1/(qtList2-qtList1)
    CASE DEFAULT
     write(*,*) 'Unexpected combination of integrations....'
    END SELECT
    
    call InitializeScaleVariations(2d0,1d0,1d0,1d0)
     SELECT CASE(Xtype)
    CASE(0)
     call CalculateXsec(xSec2,qtList1)
    CASE(1)
     call CalculateXsec_Qint(xSec2,qtList1,Q1,Q2)
    CASE(10)
     call CalculateXsec_Yint(xSec2,qtList1,y1,y2)
    CASE(100)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(11)
     call CalculateXsec_Qint_Yint(xSec2,qtList1,y1,y2,Q1,Q2)
    CASE(101)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(110)
     call CalculateXsec_PTint_Yint(xSec2,qtList1,qtList2,y1,y2)
     xSec2=xSec2/(qtList2-qtList1)
    CASE(111)
     call CalculateXsec_PTint_Qint_Yint(xSec2,qtList1,qtList2,Q1,Q2,y1,y2)
     xSec2=xSec2/(qtList2-qtList1)
    CASE DEFAULT
     write(*,*) 'Unexpected combination of integrations....'
    END SELECT
    end if
    
        If(varyMOD==2 .or. varyMOD==5 ) then
     call InitializeScaleVariations(1d0,0.5d0,1d0,1d0)
     SELECT CASE(Xtype)
    CASE(0)
     call CalculateXsec(xSec3,qtList1)
    CASE(1)
     call CalculateXsec_Qint(xSec3,qtList1,Q1,Q2)
    CASE(10)
     call CalculateXsec_Yint(xSec3,qtList1,y1,y2)
    CASE(100)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(11)
     call CalculateXsec_Qint_Yint(xSec3,qtList1,y1,y2,Q1,Q2)
    CASE(101)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(110)
     call CalculateXsec_PTint_Yint(xSec3,qtList1,qtList2,y1,y2)
     xSec3=xSec3/(qtList2-qtList1)
    CASE(111)
     call CalculateXsec_PTint_Qint_Yint(xSec3,qtList1,qtList2,Q1,Q2,y1,y2)
     xSec3=xSec3/(qtList2-qtList1)
    CASE DEFAULT
     write(*,*) 'Unexpected combination of integrations....'
    END SELECT
    
    call InitializeScaleVariations(1d0,2d0,1d0,1d0)
     SELECT CASE(Xtype)
    CASE(0)
     call CalculateXsec(xSec4,qtList1)
    CASE(1)
     call CalculateXsec_Qint(xSec4,qtList1,Q1,Q2)
    CASE(10)
     call CalculateXsec_Yint(xSec4,qtList1,y1,y2)
    CASE(100)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(11)
     call CalculateXsec_Qint_Yint(xSec4,qtList1,y1,y2,Q1,Q2)
    CASE(101)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(110)
     call CalculateXsec_PTint_Yint(xSec4,qtList1,qtList2,y1,y2)
     xSec4=xSec4/(qtList2-qtList1)
    CASE(111)
     call CalculateXsec_PTint_Qint_Yint(xSec4,qtList1,qtList2,Q1,Q2,y1,y2)
     xSec4=xSec4/(qtList2-qtList1)
    CASE DEFAULT
     write(*,*) 'Unexpected combination of integrations....'
    END SELECT
    end if
    
    
    If(varyMOD==3 .or. varyMOD==5 ) then
     call InitializeScaleVariations(1d0,1d0,0.5d0,1d0)
     SELECT CASE(Xtype)
    CASE(0)
     call CalculateXsec(xSec5,qtList1)
    CASE(1)
     call CalculateXsec_Qint(xSec5,qtList1,Q1,Q2)
    CASE(10)
     call CalculateXsec_Yint(xSec5,qtList1,y1,y2)
    CASE(100)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(11)
     call CalculateXsec_Qint_Yint(xSec5,qtList1,y1,y2,Q1,Q2)
    CASE(101)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(110)
     call CalculateXsec_PTint_Yint(xSec5,qtList1,qtList2,y1,y2)
     xSec5=xSec5/(qtList2-qtList1)
    CASE(111)
     call CalculateXsec_PTint_Qint_Yint(xSec5,qtList1,qtList2,Q1,Q2,y1,y2)
     xSec5=xSec5/(qtList2-qtList1)
    CASE DEFAULT
     write(*,*) 'Unexpected combination of integrations....'
    END SELECT
    
    call InitializeScaleVariations(1d0,1d0,2d0,1d0)
     SELECT CASE(Xtype)
    CASE(0)
     call CalculateXsec(xSec6,qtList1)
    CASE(1)
     call CalculateXsec_Qint(xSec6,qtList1,Q1,Q2)
    CASE(10)
     call CalculateXsec_Yint(xSec6,qtList1,y1,y2)
    CASE(100)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(11)
     call CalculateXsec_Qint_Yint(xSec6,qtList1,y1,y2,Q1,Q2)
    CASE(101)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(110)
     call CalculateXsec_PTint_Yint(xSec6,qtList1,qtList2,y1,y2)
     xSec6=xSec6/(qtList2-qtList1)
    CASE(111)
     call CalculateXsec_PTint_Qint_Yint(xSec6,qtList1,qtList2,Q1,Q2,y1,y2)
     xSec6=xSec6/(qtList2-qtList1)
    CASE DEFAULT
     write(*,*) 'Unexpected combination of integrations....'
    END SELECT
    end if
    
    If(varyMOD==4 .or. varyMOD==5 ) then
     call InitializeScaleVariations(1d0,1d0,1d0,0.5d0)
     SELECT CASE(Xtype)
    CASE(0)
     call CalculateXsec(xSec7,qtList1)
    CASE(1)
     call CalculateXsec_Qint(xSec7,qtList1,Q1,Q2)
    CASE(10)
     call CalculateXsec_Yint(xSec7,qtList1,y1,y2)
    CASE(100)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(11)
     call CalculateXsec_Qint_Yint(xSec7,qtList1,y1,y2,Q1,Q2)
    CASE(101)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(110)
     call CalculateXsec_PTint_Yint(xSec7,qtList1,qtList2,y1,y2)
     xSec7=xSec7/(qtList2-qtList1)
    CASE(111)
     call CalculateXsec_PTint_Qint_Yint(xSec7,qtList1,qtList2,Q1,Q2,y1,y2)
     xSec7=xSec7/(qtList2-qtList1)
    CASE DEFAULT
     write(*,*) 'Unexpected combination of integrations....'
    END SELECT
    
    call InitializeScaleVariations(1d0,1d0,1d0,2d0)
     SELECT CASE(Xtype)
    CASE(0)
     call CalculateXsec(xSec8,qtList1)
    CASE(1)
     call CalculateXsec_Qint(xSec8,qtList1,Q1,Q2)
    CASE(10)
     call CalculateXsec_Yint(xSec8,qtList1,y1,y2)
    CASE(100)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(11)
     call CalculateXsec_Qint_Yint(xSec8,qtList1,y1,y2,Q1,Q2)
    CASE(101)
     write(*,*) 'I did not expect that someone would need it. Contact me, and I will make this option.'
    CASE(110)
     call CalculateXsec_PTint_Yint(xSec8,qtList1,qtList2,y1,y2)
     xSec8=xSec8/(qtList2-qtList1)
    CASE(111)
     call CalculateXsec_PTint_Qint_Yint(xSec8,qtList1,qtList2,Q1,Q2,y1,y2)
     xSec8=xSec8/(qtList2-qtList1)
    CASE DEFAULT
     write(*,*) 'Unexpected combination of integrations....'
    END SELECT
    end if
    
   inquire(file="bin/output.dat", exist=exist)
    if (exist) then
    open(7, file="bin/output.dat", status="old", position="append", action="write")
    else
    open(7, file="bin/output.dat", status="new", action="write")
    end if
  if(varyMOD==0) then
  write(7,*) ' qT	dsigma'
  else
  write(7,*) ' qT	dsigma	-dSigma	+dSigma '
  end if
  do j=1,length
   if(qqq==0) then
    dummyREAL=qtList1(j)
   else
    dummyREAL=(qtList1(j)+qtList2(j))/2d0
   end if
  SELECT CASE(varyMOD)
    CASE(0)
     write(7,*) dummyREAL, xSec(j)
    CASE(1)
     err1=MAX(xSec(j),xSec1(j),xSec2(j))
     err2=MIN(xSec(j),xSec1(j),xSec2(j))
     write(7,*) dummyREAL, xSec(j), err2-xSec(j),err1-xSec(j)
    CASE(2)
     err1=MAX(xSec(j),xSec3(j),xSec4(j))
     err2=MIN(xSec(j),xSec3(j),xSec4(j))
     write(7,*) dummyREAL, xSec(j), err2-xSec(j),err1-xSec(j)
    CASE(3)
     err1=MAX(xSec(j),xSec5(j),xSec6(j))
     err2=MIN(xSec(j),xSec5(j),xSec6(j))
     write(7,*) dummyREAL, xSec(j), err2-xSec(j),err1-xSec(j)
    CASE(4)
     err1=MAX(xSec(j),xSec7(j),xSec8(j))
     err2=MIN(xSec(j),xSec7(j),xSec8(j))
     write(7,*) dummyREAL, xSec(j), err2-xSec(j),err1-xSec(j)
    CASE(5)
     err1=MAX(xSec(j),xSec1(j),xSec2(j),xSec3(j),xSec4(j),xSec5(j),xSec6(j),xSec7(j),xSec8(j))
     err2=MIN(xSec(j),xSec1(j),xSec2(j),xSec3(j),xSec4(j),xSec5(j),xSec6(j),xSec7(j),xSec8(j))
     write(7,*) dummyREAL, xSec(j), err2-xSec(j),err1-xSec(j)
    END SELECT
  end do
  write(7,*) ' '
  write(7,*) ' '
  write(7,*) ' '
  CLOSE (7, STATUS='KEEP')
  write(*,*) 'Evaluation finished. The result saved to output.dat'
call cpu_time(time2)

write(*,*) 'Calculation time=',time2-time1

end program DY
