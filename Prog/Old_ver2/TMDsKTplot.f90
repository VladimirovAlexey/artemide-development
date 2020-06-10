program example
use aTMDe_control
use TMDs_inKT
use TMDs
implicit none

real*8::x,mu,zeta,qT
integer::i,j,k
real*8::fTable(1:100,-5:5),fCentral(-5:5),fErr(-5:5)
character(len=100):: fname

!call artemide_Initialize('constants-file','Models/BSV19.bFIT/')
call artemide_Initialize('const-DYfit18_LO-test',prefix='/misc/data2/braun/vla18041/arTeMiDe_Repository/')
call artemide_SetReplica_uTMDPDF(0)
call artemide_SetReplica_TMDR(0)
 
 write(*,*) uTMDPDF_kT_5(0.1d0,1d0,1)
 
 do k=1,9
 
 !!! fixed x
 select case(k)
  case(1)
    x=0.01d0
    mu=5d0
   case(2)
    x=0.001d0
    mu=5d0
   case(3)
    x=0.0001d0
    mu=5d0
   case(4)
    x=0.01d0
    mu=10d0
   case(5)
    x=0.001d0
    mu=10d0
   case(6)
    x=0.0001d0
    mu=10d0
   case(7)
    x=0.01d0
    mu=100d0
   case(8)
    x=0.001d0
    mu=100d0
   case(9)
    x=0.0001d0
    mu=100d0
   case default
    write(*,*) 'YOU IDIOT!!!'
    stop
   end select
 zeta=mu**2
 
 write(*,*) 'Case ',x,mu
 
 write(fname,"('fX_',I1,'_KT.dat')") k
 OPEN(UNIT=52, FILE=trim(fname), ACTION="write", STATUS="replace")
 
 write(52,*) 'unpolarized TMDPDF in KT,		artemide v2.00'
 write(52,*) 'fit =    BSV19.bFIT'
 write(52,"('x=',F6.4)"),x
 write(52,"('mu=',F4.0)"),mu
 write(52,*) 'qT, bBar, bBar-, bBar+, cBar, cBar-, cBar+, sBar, sBar-, sBar+, uBar, uBar-, uBar+, dBar, dBar-, dBar+,&
	      gluon=0,gluon-=0,gluon+=0, d, d-, d+, u ,u-, u+, s, s-, s+, c, c-, c+, b, b-, b+'
 
 do i=0,30
  qT=0.1d0-1d0+10**(0.05d0*i)
  do j=1,100
    call artemide_SetReplica_uTMDPDF(j)
    call artemide_SetReplica_TMDR(j)
    fTable(j,-5:5)=x*uTMDPDF_kT_5(x,qT,mu,zeta,1)
  end do
  
  fCentral=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  do j=1,100
    fCentral=fCentral+fTable(j,-5:5)
  end do
  fCentral=fCentral/100d0
  
  fErr=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  do j=1,100
    fErr=fErr+(fTable(j,-5:5)-fCentral)**2
  end do
  fErr=sqrt(fErr/100d0)
  
  
  write(52,"('  ',F10.4, '  ')",advance='no') qT
  do j=-5,5
    write(52,"('  ',ES12.4, '  ')",advance='no') fCentral(j)
    write(52,"('  ',ES12.4, '  ')",advance='no') fCentral(j)-fErr(j)
    write(52,"('  ',ES12.4, '  ')",advance='no') fCentral(j)+fErr(j)
  end do
  write(52,*) ' '
  write(*,*) i
 end do
 CLOSE (52, STATUS='KEEP')
 end do
 
 
 !---------------------------------------------------------------------------------------------------------------------------
 
 write(*,*) '=================================================================='
 write(*,*) '=================================================================='
 write(*,*) '=================================================================='
 
 
 do k=1,9
 
 !!! fixed x
 select case(k)
  case(1)
    qT=3d0
    mu=5d0
   case(2)
    qT=5d0
    mu=5d0
   case(3)
    qT=10d0
    mu=5d0
   case(4)
    qT=3d0
    mu=10d0
   case(5)
    qT=5d0
    mu=10d0
   case(6)
    qT=10d0
    mu=10d0
   case(7)
    qT=3d0
    mu=100d0
   case(8)
    qT=5d0
    mu=100d0
   case(9)
    qT=10d0
    mu=100d0
   case default
    write(*,*) 'YOU IDIOT!!!'
    stop
   end select
 zeta=mu**2
 
 write(*,*) 'Case ',x,mu
 
 write(fname,"('fKT_',I1,'_X.dat')") k
 OPEN(UNIT=52, FILE=trim(fname), ACTION="write", STATUS="replace")
 
 write(52,*) 'unpolarized TMDPDF in KT,		artemide v2.00'
 write(52,*) 'fit =    BSV19.bFIT'
 write(52,"('qT=',F6.1)"),qT
 write(52,"('mu=',F4.0)"),mu
 write(52,*) 'x, bBar, bBar-, bBar+, cBar, cBar-, cBar+, sBar, sBar-, sBar+, uBar, uBar-, uBar+, dBar, dBar-, dBar+,&
	      gluon=0,gluon-=0,gluon+=0, d, d-, d+, u ,u-, u+, s, s-, s+, c, c-, c+, b, b-, b+'
 
 do i=0,50
  x=10**(-0.08d0*i)
  do j=1,100
    call artemide_SetReplica_uTMDPDF(j)
    call artemide_SetReplica_TMDR(j)
    fTable(j,-5:5)=x*uTMDPDF_kT_5(x,qT,mu,zeta,1)
  end do
  
  fCentral=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  do j=1,100
    fCentral=fCentral+fTable(j,-5:5)
  end do
  fCentral=fCentral/100d0
  
  fErr=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  do j=1,100
    fErr=fErr+(fTable(j,-5:5)-fCentral)**2
  end do
  fErr=sqrt(fErr/100d0)
  
  
  write(52,"('  ',F10.4, '  ')",advance='no') x
  do j=-5,5
    write(52,"('  ',ES12.4, '  ')",advance='no') fCentral(j)
    write(52,"('  ',ES12.4, '  ')",advance='no') fCentral(j)-fErr(j)
    write(52,"('  ',ES12.4, '  ')",advance='no') fCentral(j)+fErr(j)
  end do
  write(52,*) ' '
  write(*,*) i
 end do
 CLOSE (52, STATUS='KEEP')
 end do

end program example