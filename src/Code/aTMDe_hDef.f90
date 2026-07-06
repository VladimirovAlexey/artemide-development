!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       arTeMiDe 3.05
!
!   Contains definition of transformation of parton content of hadrons
!
!                               A.Vladimirov (06.07.2026)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module aTMDe_hDef
use aTMDe_numerics
use aTMDe_IO
implicit none
private

public::PDFforH!,transformFFtoh

contains

!!!!! Takes input pdf vector and transforms it for the target hadron
pure function PDFforH(pdf,h)
real(dp),dimension(-5:5),intent(in)::pdf
integer,intent(in)::h
real(dp),dimension(-5:5)::PDFforH

integer::A,Z
real(dp)::za,ud,ud_bar

SELECT CASE(abs(h))
    CASE(11)
        !!!!! neutron case (u<->d)
        PDFforH=[pdf(-5),pdf(-4),pdf(-3),pdf(-1),pdf(-2),pdf(0),pdf(2),pdf(1),pdf(3),pdf(4),pdf(5)]
    CASE(12)
        !!!!! deuteron case (p+n)/2
        ud_bar=0.5*(pdf(-1)+pdf(-2))
        ud=0.5*(pdf(1)+pdf(2))
        PDFforH=[pdf(-5),pdf(-4),pdf(-3),ud_bar,ud_bar,&
                pdf(0),ud,ud,pdf(3),pdf(4),pdf(5)]
    CASE(64029)
        !!! copper A=64, Z=29, N=35
        PDFforH=[pdf(-5),pdf(-4),pdf(-3),29._dp/64._dp*pdf(-2)+35._dp/64._dp*pdf(-1),29._dp/64._dp*pdf(-1)+35._dp/64._dp*pdf(-2),&
                pdf(0),&
                29._dp/64._dp*pdf(1)+35._dp/64._dp*pdf(2),29._dp/64._dp*pdf(2)+35._dp/64._dp*pdf(1),pdf(3),pdf(4),pdf(5)]
    CASE(184074)
        !!! tungsten A=184, Z=74, N=110
        PDFforH=[pdf(-5),pdf(-4),pdf(-3),37._dp/92._dp*pdf(-2)+55._dp/92._dp*pdf(-1),37._dp/92._dp*pdf(-1)+55._dp/92._dp*pdf(-2),&
                pdf(0),&
                37._dp/92._dp*pdf(1)+55._dp/92._dp*pdf(2),37._dp/92._dp*pdf(2)+55._dp/92._dp*pdf(1),pdf(3),pdf(4),pdf(5)]

    CASE DEFAULT
        !!! very unusual cases (generally these cases should not appear)
        if(abs(h)<10) then
            !!!!! direct pass (it should be checked before entry to this routine)
            PDFforH(-5:5)=pdf(-5:5)
        else if(abs(h)>1000 .and. abs(h)<1000000) then !!! general nuclei case
            !!!!! for frequently used nuclei, better to hard-code them as specific CASE branches above
            A=abs(h)/1000
            Z=mod(abs(h),1000)
            za=real(Z,dp)/real(A,dp)
            PDFforH=[pdf(-5),pdf(-4),pdf(-3),za*pdf(-2)+(1._dp-za)*pdf(-1),za*pdf(-1)+(1._dp-za)*pdf(-2),&
                pdf(0),&
                za*pdf(1)+(1._dp-za)*pdf(2),za*pdf(2)+(1._dp-za)*pdf(1),pdf(3),pdf(4),pdf(5)]
        else
            error stop ErrorString("Unknown type of hadron for PDF "//trim(numToStr(h)),"aTMDe-hDef")
        end if
END SELECT

!!!! finally revert the anti-hadron content
if(h<0) PDFforH=PDFforH(5:-5:-1)

end function PDFforH

end module aTMDe_hDef
