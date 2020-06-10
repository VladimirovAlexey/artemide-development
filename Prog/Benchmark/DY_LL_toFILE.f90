!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDX_DY
implicit none

character(*),parameter::location='/misc/data2/braun/vla18041/WorkingFiles/TMD/arTeMiDe/Benchmark/aTMDE-tables-level1/'
character(*),parameter::order='N3LL'
character(4)::Qname

integer,dimension(1:3)::proc
real*8::y,Q,s
real*8::stepqT,maxqT,r
integer::i,maxI
real*8,allocatable::X(:),qT(:),X1(:),X2(:),X3(:),X4(:)

call artemide_Initialize('const-Bench_'//trim(order),&
	prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/Constants-files/DY_Benchmark/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call artemide_SetNPparameters_TMDR((/2.86041d0, 0.229551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.05d0, 0.05d0, 0d0, 0d0, 0d0, 0.d0, 0.d0/))

! call artemide_SetNPparameters_TMDR((/2d0,0.001d0/))
! call artemide_SetNPparameters_uTMDPDF((/0.01d0, 0.01d0, 0d0, 0d0, 0d0, 0.d0, 0.d0/))


proc=(/1,1,5/) !!!! pp DY
s=13000d0**2
Q=91.15348061918276d0

maxI=120!Int(maxqT/stepqT)
allocate(qT(1:maxI))
allocate(X(1:maxI))

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
    write(*,*) '------------------------------Y=0 -------------------------------'
    y=0.0d0
    call TMDX_DY_XSetup(s,Q,y)

    call CalcXsec_DY(X,qT)

    OPEN(UNIT=55, FILE=location//"aTMDe_"//trim(Qname)//"_Y0_"//trim(order)//"_level1.txt", ACTION="write", STATUS="replace")
    write(55,*) "# artemide results for canonical resummation (CSS-type1)"
    write(55,*) "# s=13000 GeV"
    write(55,*) "# Q=",Q," GeV"
    write(55,*) "# Y=",y
    write(55,*) "# order=",order
    write(55,*) "#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY"
    do i=1,maxI
        write(55,'(F6.2,"        ",F14.12)') qT(i),2d0*Q*2d0*qT(i)*X(i)
    end do
    CLOSE(55,STATUS='KEEP')

    write(*,*) '------------------------------Y=1.2-------------------------------'
    y=1.2d0
    call TMDX_DY_XSetup(s,Q,y)

    call CalcXsec_DY(X,qT)

    OPEN(UNIT=55, FILE=location//"aTMDe_"//trim(Qname)//"_Y1.2_"//trim(order)//"_level1.txt", ACTION="write", STATUS="replace")
    write(55,*) "# artemide results for canonical resummation (CSS-type1)"
    write(55,*) "# s=13000 GeV"
    write(55,*) "# Q=",Q," GeV"
    write(55,*) "# Y=",y
    write(55,*) "# order=",order
    write(55,*) "#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY"
    do i=1,maxI
        write(55,'(F6.2,"        ",F14.12)') qT(i),2d0*Q*2d0*qT(i)*X(i)
    end do
    CLOSE(55,STATUS='KEEP')

    write(*,*) '------------------------------Y=2.4-------------------------------'
    y=2.4d0
    call TMDX_DY_XSetup(s,Q,y)

    call CalcXsec_DY(X,qT)

    OPEN(UNIT=55, FILE=location//"aTMDe_"//trim(Qname)//"_Y2.4_"//trim(order)//"_level1.txt", ACTION="write", STATUS="replace")
    write(55,*) "# artemide results for canonical resummation (CSS-type1)"
    write(55,*) "# s=13000 GeV"
    write(55,*) "# Q=",Q," GeV"
    write(55,*) "# Y=",y
    write(55,*) "# order=",order
    write(55,*) "#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY"
    do i=1,maxI
        write(55,'(F6.2,"        ",F14.12)') qT(i),2d0*Q*2d0*qT(i)*X(i)
    end do
    CLOSE(55,STATUS='KEEP')

    write(*,*) '------------------------------Y=3.6-------------------------------'
    y=3.6d0
    call TMDX_DY_XSetup(s,Q,y)

    call CalcXsec_DY(X,qT)

    OPEN(UNIT=55, FILE=location//"aTMDe_"//trim(Qname)//"_Y3.6_"//trim(order)//"_level1.txt", ACTION="write", STATUS="replace")
    write(55,*) "# artemide results for canonical resummation (CSS-type1)"
    write(55,*) "# s=13000 GeV"
    write(55,*) "# Q=",Q," GeV"
    write(55,*) "# Y=",y
    write(55,*) "# order=",order
    write(55,*) "#qT   |     2qT 2Q dsigma/dqT^2 /dQ^2/dY"
    do i=1,maxI
        write(55,'(F6.2,"        ",F14.12)') qT(i),2d0*Q*2d0*qT(i)*X(i)
    end do
    CLOSE(55,STATUS='KEEP')
end subroutine SaveDifferentY

end program example