!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes reads the file with (number ,chi^2) and plot PDF-replicas as function of x.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use QCDinput
implicit none

integer::numR=1000
integer::i,j,k
real*8::chi2,inputPDF(-6:6)
real*8::x

call InitPDFsetByName("NNPDF31_nnlo_as_0118_1000")
call InitPDF(0)

call evolvePDF(0.1d0,2d0,inputPDF)
      
write(*,*) inputPDF(-2:2)


OPEN(UNIT=51, FILE="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/FIGURES_DY+SIDIS_2019/PDFrep/Chi2_for_PDFrep",&
      ACTION="read", STATUS="old")
 !! replicas
 !numR=100
OPEN(UNIT=52, FILE="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/FIGURES_DY+SIDIS_2019/PDFrep/PDFplot_rep_d.dat",&
      ACTION="write", STATUS="new")
OPEN(UNIT=53, FILE="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/FIGURES_DY+SIDIS_2019/PDFrep/PDFplot_rep_u.dat",&
      ACTION="write", STATUS="new")
OPEN(UNIT=54, FILE="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/FIGURES_DY+SIDIS_2019/PDFrep/PDFplot_rep_s.dat",&
      ACTION="write", STATUS="new")
OPEN(UNIT=55, FILE="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/FIGURES_DY+SIDIS_2019/PDFrep/PDFplot_rep_g.dat",&
      ACTION="write", STATUS="new")

write(52,'("{")',advance="no") 
write(53,'("{")',advance="no") 
write(54,'("{")',advance="no") 
write(55,'("{")',advance="no") 
do j=1,numR+1
    read(51,*) k,chi2
    call InitPDF(k)
    write(52,'("{",F16.10,",{")',advance="no") chi2
    write(53,'("{",F16.10,",{")',advance="no") chi2
    write(54,'("{",F16.10,",{")',advance="no") chi2
    write(55,'("{",F16.10,",{")',advance="no") chi2
    do i=0,99
        x=10d0**(-4d0+4d0*i/100d0)
        call evolvePDF(x,2d0,inputPDF)  
        write(52,'("{",F8.5,",",F16.10,"},")',advance="no") x,inputPDF(1)
        write(53,'("{",F8.5,",",F16.10,"},")',advance="no") x,inputPDF(2)
        write(54,'("{",F8.5,",",F16.10,"},")',advance="no") x,inputPDF(3)
        write(55,'("{",F8.5,",",F16.10,"},")',advance="no") x,inputPDF(0)
    end do
    write(52,'("{1., 0.}}},")',advance="no")
    write(53,'("{1., 0.}}},")',advance="no")   
    write(54,'("{1., 0.}}},")',advance="no")   
    write(55,'("{1., 0.}}},")',advance="no")   
end do
write(52,'("{0., 0.}}")',advance="no")
write(53,'("{0., 0.}}")',advance="no")   
write(54,'("{0., 0.}}")',advance="no")   
write(55,'("{0., 0.}}")',advance="no")   
CLOSE (51, STATUS='KEEP')
CLOSE (52, STATUS='KEEP')
CLOSE (53, STATUS='KEEP')
CLOSE (54, STATUS='KEEP')
CLOSE (55, STATUS='KEEP')
end program example
