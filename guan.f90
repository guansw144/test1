module guan
use laser
implicit none

contains

subroutine FFW(nn,var1,var2,sigd)
implicit none
integer::nn,sigd
double complex::var1(nn),var2(nn)
integer*8::plan

call dfftw_plan_dft_1d(plan,nn,var1,var2,sigd,FFTW_ESTIMATE)
call dfftw_execute_dft(plan,var1,var2)
call dfftw_destroy_plan(plan)

end subroutine
!----------------------------------------------------------------------------------------------------
real(8) function eev(kxx,kyy)
implicit none
real(8)::kxx,kyy
eev=-vf*(kxx*kxx+kyy*kyy)**0.5d0
end function
!---------------------------------------------------------------------------------------------------------
real(8) function eec(kxx,kyy)
implicit none
real(8)::kxx,kyy

eec=vf*(kxx*kxx+kyy*kyy)**0.5d0
end function
!--------------------------------------------------------------------------------------------------
real(8) function ddx(kxx,kyy)
implicit none
real(8)::kxx,kyy

ddx=0.5d0*kyy/(kxx*kxx+kyy*kyy)
end function
!-----------------------------------------------------------------------------------------------------
real(8) function ddy(kxx,kyy)
implicit none
real(8)::kxx,kyy

ddy=-0.5d0*kxx/(kxx*kxx+kyy*kyy)
end function
!---------------------------------------------------------------------------------------------------------
real(8) function themm(kxx,kyy)
implicit none
real(8)::kxx,kyy

themm=atan2(kyy,kxx)
end function
!------------------------------------------------------------------
subroutine f02 ( t, y, yp )
implicit none
real(8)::t
real(8)::y(nvar)
real(8)::yp(nvar)
common px,py
real(8)::px,py
double complex::sdx,sdy,dpps,ap1,ap2,ap3
real(8)::egv,egc

sdx=ddx(px+aax(t),py+aay(t))
sdy=ddy(px+aax(t),py+aay(t))
dpps=2.d0*(eex(t)*sdx+eey(t)*sdy)
egv=eev(px+aax(t),py+aay(t))
egc=eec(px+aax(t),py+aay(t))

ap1=y(1)+ui*y(2)
ap2=y(3)+ui*y(4)
ap3=y(5)+ui*y(6)

yp(1)=dreal(-ui*(egc-egv)*ap1+ui*eex(t)*sdx*(1-ap2-ap3)+ui*eey(t)*sdy*(1-ap2-ap3)-tss1*ap1)
yp(2)=dimag(-ui*(egc-egv)*ap1+ui*eex(t)*sdx*(1-ap2-ap3)+ui*eey(t)*sdy*(1-ap2-ap3)-tss1*ap1)
yp(3)=dreal(ui*eex(t)*(sdx*conjg(ap1)-conjg(sdx)*ap1)+ui*eey(t)*(sdy*conjg(ap1)-conjg(sdy)*ap1))
yp(4)=dimag(ui*eex(t)*(sdx*conjg(ap1)-conjg(sdx)*ap1)+ui*eey(t)*(sdy*conjg(ap1)-conjg(sdy)*ap1))
yp(5)=dreal(ui*eex(t)*(sdx*conjg(ap1)-conjg(sdx)*ap1)+ui*eey(t)*(sdy*conjg(ap1)-conjg(sdy)*ap1))
yp(6)=dimag(ui*eex(t)*(sdx*conjg(ap1)-conjg(sdx)*ap1)+ui*eey(t)*(sdy*conjg(ap1)-conjg(sdy)*ap1))

  return
end  subroutine
end module