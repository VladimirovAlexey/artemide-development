!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDX_DY
implicit none

character(*),parameter::location='/home/vla18041/LinkData2/WorkingFiles/TMD/arTeMiDe/Benchmark/aTMDe-tables-level2/'
!character(*),parameter::location='/home/alexey/WorkingFiles/TMD/Benchmark/Level2/'
character(*),parameter::order='NLL'
character(4)::Qname

integer,dimension(1:3)::proc
real*8::y,Q,s
real*8::stepqT,maxqT,r
integer::i,maxI
real*8,allocatable::qT(:),Xa(:),Xa1(:),Xa2(:),Xa3(:),Xa4(:),&
                            Xb(:),Xb1(:),Xb2(:),Xb3(:),Xb4(:),&
                            Xc(:),Xc1(:),Xc2(:),Xc3(:),Xc4(:),&
                            Xd(:),Xd1(:),Xd2(:),Xd3(:),Xd4(:)

call artemide_Initialize('const-Bench_'//trim(order)//'-type3',&
	prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/Constants-files/DY_Benchmark/')
	!prefix='/home/alexey/artemide_Repository/Constants-files/DY_Benchmark/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call artemide_SetNPparameters_TMDR((/1.86d0, 0.0296d0/))
call artemide_SetNPparameters_uTMDPDF((/0.005d0, 0.005d0, 0d0, 0d0, 0d0, 0.d0, 0.d0/))

! call artemide_SetNPparameters_TMDR((/2d0,0.001d0/))
! call artemide_SetNPparameters_uTMDPDF((/0.01d0, 0.01d0, 0d0, 0d0, 0d0, 0.d0, 0.d0/))


proc=(/1,1,5/) !!!! pp DY
s=13000d0**2
Q=91.15348061918276d0

maxI=120!Int(maxqT/stepqT)
allocate(qT(1:maxI))
allocate(Xa(1:maxI),Xa1(1:maxI),Xa2(1:maxI),Xa3(1:maxI),Xa4(1:maxI))
allocate(Xb(1:maxI),Xb1(1:maxI),Xb2(1:maxI),Xb3(1:maxI),Xb4(1:maxI))
allocate(Xc(1:maxI),Xc1(1:maxI),Xc2(1:maxI),Xc3(1:maxI),Xc4(1:maxI))
allocate(Xd(1:maxI),Xd1(1:maxI),Xd2(1:maxI),Xd3(1:maxI),Xd4(1:maxI))


!--------------------------------qT-table
do i=1,40
  qT(i)=0.25d0*i
end do
do i=41,80
  qT(i)=0.5d0*(i-40)+10
end do
do i=81,100
  qT(i)=1d0*(i-80)+30
end do
do i=101,maxI
  qT(i)=2.5d0*(i-100)+50
end do


call TMDX_DY_setProcess(proc)
call TMDX_DY_XSetup(s,Q,y)
call TMDX_DY_SetCuts(.false.,0d0,-10d0,10d0)!!! no cuts



Qname="QmZ"
Q=91.15348061918276d0
write(*,*)'>>>>>>>>>>>>>>>>>>>>',Qname,'.....',Q,'<<<<<<<<<<<<<<<<<<<<'
call SaveDifferentY()

Qname="66"
Q=66d0
write(*,*)'>>>>>>>>>>>>>>>>>>>>',Qname,'.....',Q,'<<<<<<<<<<<<<<<<<<<<'
call SaveDifferentY()

Qname="116"
Q=116d0
write(*,*)'>>>>>>>>>>>>>>>>>>>>',Qname,'.....',Q,'<<<<<<<<<<<<<<<<<<<<'
call SaveDifferentY()

Qname="300"
Q=300d0
write(*,*)'>>>>>>>>>>>>>>>>>>>>',Qname,'.....',Q,'<<<<<<<<<<<<<<<<<<<<'
call SaveDifferentY()

Qname="1000"
Q=1000d0
write(*,*)'>>>>>>>>>>>>>>>>>>>>',Qname,'.....',Q,'<<<<<<<<<<<<<<<<<<<<'
call SaveDifferentY()

contains 

subroutine SaveDifferentY()
    real*8, parameter:: r=2d0
    
    write(*,*) '------------------------------CENTRAL VALUE-------------------------------'    
    call artemide_SetScaleVariations(1d0,1d0,1d0,1d0)
    
    y=0.0d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xa,qT)
    
    y=1.2d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xb,qT)
    
    y=2.4d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xc,qT)
    
    y=3.6d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xd,qT)
    
    write(*,*) '------------------------------c2- VALUE-------------------------------'    
    call artemide_SetScaleVariations(1d0,1d0/r,1d0,1d0)
    
    y=0.0d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xa1,qT)
    
    y=1.2d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xb1,qT)
    
    y=2.4d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xc1,qT)
    
    y=3.6d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xd1,qT)
    
    write(*,*) '------------------------------c2+ VALUE-------------------------------'    
    call artemide_SetScaleVariations(1d0,1d0*r,1d0,1d0)
    
    y=0.0d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xa2,qT)
    
    y=1.2d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xb2,qT)
    
    y=2.4d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xc2,qT)
    
    y=3.6d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xd2,qT)
    
    write(*,*) '------------------------------c4- VALUE-------------------------------'    
    call artemide_SetScaleVariations(1d0,1d0,1d0,1d0/r)
    
    y=0.0d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xa3,qT)
    
    y=1.2d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xb3,qT)
    
    y=2.4d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xc3,qT)
    
    y=3.6d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xd3,qT)
    
    write(*,*) '------------------------------c4+ VALUE-------------------------------'    
    call artemide_SetScaleVariations(1d0,1d0,1d0,1d0*r)
    
    y=0.0d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xa4,qT)
    
    y=1.2d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xb4,qT)
    
    y=2.4d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xc4,qT)
    
    y=3.6d0
    call TMDX_DY_XSetup(s,Q,y)
    call CalcXsec_DY(Xd4,qT)
  
    write(*,*) '------------------------------SAVING RESULT-------------------------------'    


    OPEN(UNIT=55, FILE=location//"aTMDe_"//trim(Qname)//"_Y0_"//trim(order)//"_level2.txt", ACTION="write", STATUS="replace")
    !write(55,*) "# artemide results for canonical resummation (CSS-type1)"
    write(55,*) "# artemide results for level-2 comparison (default evolution settings)"
    write(55,*) "# s=13000 GeV"
    write(55,*) "# Q=",Q," GeV"
    write(55,*) "# Y=",y
    write(55,*) "# order=",order
    write(55,*) "#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY     |  -scale var.  |   +scale var."
    do i=1,maxI        
        write(55,'(F6.2,"        ",F18.12,"        ",F18.12,"        ",F18.12)') qT(i),2d0*Q*2d0*qT(i)*Xa(i), &
                2d0*Q*2d0*qT(i)*(min(Xa(i),Xa1(i),Xa2(i),Xa3(i),Xa4(i))-Xa(i)), &
                2d0*Q*2d0*qT(i)*(max(Xa(i),Xa1(i),Xa2(i),Xa3(i),Xa4(i))-Xa(i))
    end do
    CLOSE(55,STATUS='KEEP')

    
    OPEN(UNIT=55, FILE=location//"aTMDe_"//trim(Qname)//"_Y1.2_"//trim(order)//"_level2.txt", ACTION="write", STATUS="replace")
    !write(55,*) "# artemide results for canonical resummation (CSS-type1)"
    write(55,*) "# artemide results for level-2 comparison (default evolution settings)"
    write(55,*) "# s=13000 GeV"
    write(55,*) "# Q=",Q," GeV"
    write(55,*) "# Y=",y
    write(55,*) "# order=",order
    write(55,*) "#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY     |  -scale var.  |   +scale var."
    do i=1,maxI        
        write(55,'(F6.2,"        ",F18.12,"        ",F18.12,"        ",F18.12)') qT(i),2d0*Q*2d0*qT(i)*Xb(i), &
                2d0*Q*2d0*qT(i)*(min(Xb(i),Xb1(i),Xb2(i),Xb3(i),Xb4(i))-Xb(i)), &
                2d0*Q*2d0*qT(i)*(max(Xb(i),Xb1(i),Xb2(i),Xb3(i),Xb4(i))-Xb(i))
    end do
    CLOSE(55,STATUS='KEEP')

    
    OPEN(UNIT=55, FILE=location//"aTMDe_"//trim(Qname)//"_Y2.4_"//trim(order)//"_level2.txt", ACTION="write", STATUS="replace")
    !write(55,*) "# artemide results for canonical resummation (CSS-type1)"
    write(55,*) "# artemide results for level-2 comparison (default evolution settings)"
    write(55,*) "# s=13000 GeV"
    write(55,*) "# Q=",Q," GeV"
    write(55,*) "# Y=",y
    write(55,*) "# order=",order
    write(55,*) "#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY     |  -scale var.  |   +scale var."
    do i=1,maxI        
        write(55,'(F6.2,"        ",F18.12,"        ",F18.12,"        ",F18.12)') qT(i),2d0*Q*2d0*qT(i)*Xc(i), &
                2d0*Q*2d0*qT(i)*(min(Xc(i),Xc1(i),Xc2(i),Xc3(i),Xc4(i))-Xc(i)), &
                2d0*Q*2d0*qT(i)*(max(Xc(i),Xc1(i),Xc2(i),Xc3(i),Xc4(i))-Xc(i))
    end do
    CLOSE(55,STATUS='KEEP')


    OPEN(UNIT=55, FILE=location//"aTMDe_"//trim(Qname)//"_Y3.6_"//trim(order)//"_level2.txt", ACTION="write", STATUS="replace")
    !write(55,*) "# artemide results for canonical resummation (CSS-type1)"
    write(55,*) "# artemide results for level-2 comparison (default evolution settings)"
    write(55,*) "# s=13000 GeV"
    write(55,*) "# Q=",Q," GeV"
    write(55,*) "# Y=",y
    write(55,*) "# order=",order
    write(55,*) "#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY     |  -scale var.  |   +scale var."
    do i=1,maxI        
        write(55,'(F6.2,"        ",F18.12,"        ",F18.12,"        ",F18.12)') qT(i),2d0*Q*2d0*qT(i)*Xd(i), &
                2d0*Q*2d0*qT(i)*(min(Xd(i),Xd1(i),Xd2(i),Xd3(i),Xd4(i))-Xd(i)), &
                2d0*Q*2d0*qT(i)*(max(Xd(i),Xd1(i),Xd2(i),Xd3(i),Xd4(i))-Xd(i))
    end do
    CLOSE(55,STATUS='KEEP')
    
end subroutine SaveDifferentY

end program example
