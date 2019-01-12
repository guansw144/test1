program main
use guan
implicit none
integer(kind=4)::i,j,ix,ixp,it,jt,iflag,ip,ikex
real(8)::abserr, relerr
integer(kind=4)::iwork(5)
real(8)::tx,th,thh,th6,t
real(8)::tout
real(8)::work(100+21*nvar)
real(8)::y(nvar)
common px,py
real(8)::px,py
real(8)::tp(ntmsax),ts(ntmax)
real(8)::flx(ntmsax),fly(ntmsax),elx(ntmsax),ely(ntmsax)
real(8)::t0,tcenter,thetm,sins,coss
double complex::ak1,ak2,ak3,nss,ros
double complex::jjx(ntmax),jjy(ntmax),sumppx(ntmax),sumppy(ntmax)
real(8)::dthk(nnk),ths(nnk,nnt*nnk)
real(8)::ppx(nnk,nnt*nnk),ppy(nnk,nnt*nnk)
real(8)::psx(nnk,nnt*nnk),psy(nnk,nnt*nnk)
real(8)::dee,pett(1:nnk),pstt(1:nnk)
real(8)::sdkx((1+nnk)*nnk*nnt/2),sdky((1+nnk)*nnk*nnt/2)
real(8)::dww,ww(ntmax)
real(8)::jxx(ntmax),jyy(ntmax),jtot(ntmax)
double complex::jpp1(ntmax),jpp2(ntmax),aclxx(ntmax),aclyy(ntmax)
integer::comm_sz,ierr,my_rank
integer::local_a,local_b,local_n
integer::source
real(8)::local_intx(ntmax),local_inty(ntmax)
real(8)::total_intx(ntmax),total_inty(ntmax)
!-------------------------------------------------

!call MPI_INIT(ierr)
!call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
!call MPI_COMM_SIZE(MPI_COMM_WORLD,comm_sz,ierr)

open(70,file='mesh.dat',status='unknown')
open(61,file='xcurrent.dat',status='unknown')
open(62,file='ycurrent.dat',status='unknown')
open(301,file='hhgx.dat',status='unknown')
open(302,file='hhgy.dat',status='unknown')
open(303,file='hhgt.dat',status='unknown')


dee=-4.4d0*vf*aa0/dble(nnk)

do ix=1,nnk
pett(ix)=dble(ix)*dee
pstt(ix)=abs(pett(ix))/vf
end do

do ix=1,nnk
dthk(ix)=2*ppi/dble(nnt)/dble(ix)
do ixp=1,nnt*ix
ths(ix,ixp)=(ixp-1)*dthk(ix)
ppx(ix,ixp)=pstt(ix)*dcos(ths(ix,ixp))
ppy(ix,ixp)=pstt(ix)*dsin(ths(ix,ixp))
end do
end do

do ix=1,nnk
do ixp=1,nnt*ix-1
psx(ix,ixp)=(ppx(ix,ixp)+ppx(ix,ixp+1))/2.d0
psy(ix,ixp)=(ppy(ix,ixp)+ppy(ix,ixp+1))/2.d0
end do
psx(ix,nnt*ix)=(ppx(ix,nnt*ix)+ppx(ix,1))/2.d0
psy(ix,nnt*ix)=(ppy(ix,nnt*ix)+ppy(ix,1))/2.d0
end do


ip=0

do ix=1,nnk
do ixp=1,nnt*ix
ip=ip+1
sdkx(ip)=psx(ix,ixp)
sdky(ip)=psy(ix,ixp)
end do
end do

do ix=1,((1+nnk)*nnk*nnt/2)
write(70,*)sdkx(ix),sdky(ix)
end do



call asx(t0,tcenter)

do i=1,ntmsax
tp(i)=t0+dble(i-1)*dts
end do

do i=1,ntmax
ts(i)=t0+dble(i-1)*dt
end do

do it=1,ntmsax
flx(it)=aax(tp(it))
fly(it)=aay(tp(it))
elx(it)=eex(tp(it))
ely(it)=eey(tp(it))
!write(11,*)tp(it),elx(it),ely(it)
!write(12,*)tp(it),flx(it),fly(it)
end do

!--------------------------------------------------------------------------------
dww=2*ppi/(dble(ntmax)*dt)

do i=1,ntmax/2+1
ww(i)=dww*dble(i-1)
end do

ikex=ntmax/2+1

do j=ntmax/2+2,ntmax
ikex=ikex-1
ww(j)=-ww(ikex)
end do

!----------------------------------------------------------------------------------
sumppx=(0.d0,0.d0)
sumppy=(0.d0,0.d0)

local_n=((1+nnk)*nnk*nnt/2)/comm_sz
local_a=my_rank*local_n+1
local_b=local_a+local_n-1

!do ix=local_a,local_b
do ix=1,((1+nnk)*nnk*nnt/2)
px=sdkx(ix)
py=sdky(ix)

t=t0

y(1)=0.d0
y(2)=0.d0
y(3)=0.d0
y(4)=0.d0
y(5)=0.d0
y(6)=0.d0

jjx=(0.d0,0.d0)
jjy=(0.d0,0.d0)


do i=1,ntmax

thetm=themm(px+aax(t),py+aay(t))

ak1=y(1)+ui*y(2)
ak2=y(3)+ui*y(4)
ak3=y(5)+ui*y(6)


nss=ak2+ak3-1
ros=ak1

sins=dsin(thetm)
coss=dcos(thetm)


jjx(i)=nss*coss+ui*sins*(ros-conjg(ros))
jjy(i)=nss*sins-ui*coss*(ros-conjg(ros))

tout=t+dt

abserr=1.d-12
relerr=1.d-12
iflag=1

call ode(f02,nvar,y,t,tout,relerr,abserr,iflag,work,iwork)

write(100,*)iflag

if(iflag /=2) then
write(*,'(a)')' '
write(*,'(a)')'test-Fatal error!'
 stop
end if
end do

do it=1,ntmax
sumppx(it)=sumppx(it)+jjx(it)
sumppy(it)=sumppy(it)+jjy(it)
end do

end do

do it=1,ntmax
    write(61,*)ts(it)/period,dreal(sumppx(it)),dimag(sumppx(it))
    write(62,*)ts(it)/period,dreal(sumppy(it)),dimag(sumppy(it))
end do


!do i=1,ntmax
!local_intx(i)=dreal(sumppx(i))
!local_inty(i)=dreal(sumppy(i))
!end do
!
!if(my_rank/=0)then
!call MPI_send(local_intx,ntmax,mpi_double,0,0,mpi_comm_world,ierr)
!call MPI_send(local_inty,ntmax,mpi_double,0,1,mpi_comm_world,ierr)
!
!else
!
!do i=1,ntmax
!total_intx(i)=local_intx(i)
!total_inty(i)=local_inty(i)
!end do
!
!do source=1,comm_sz-1
!
!call MPI_recv(local_intx,ntmax,mpi_double,source,0,mpi_comm_world,mpi_status_ignore,ierr)
!call MPI_recv(local_inty,ntmax,mpi_double,source,1,mpi_comm_world,mpi_status_ignore,ierr)
!
!do i=1,ntmax
!total_intx(i)=total_intx(i)+local_intx(i)
!total_inty(i)=total_inty(i)+local_inty(i)
!end do
!end do
!end if
!
!if(my_rank==0)then
!
!do i=1,ntmax
!write(61,*)ts(i)/period,total_intx(i)/((1+nnk)*nnk*nnt/2)
!write(62,*)ts(i)/period,total_inty(i)/((1+nnk)*nnk*nnt/2)
!end do
!
!aclxx=total_intx*ur/((1+nnk)*nnk*nnt/2)
!aclyy=total_inty*ur/((1+nnk)*nnk*nnt/2)
!
!call FFW(ntmax,aclxx,jpp1,-1)
!call FFW(ntmax,aclyy,jpp2,-1)
!
!do it=1,ntmax
!jxx(it)=(dreal(jpp1(it))**2+dimag(jpp1(it))**2)/dble(ntmax)
!jyy(it)=(dreal(jpp2(it))**2+dimag(jpp2(it))**2)/dble(ntmax)
!jtot(it)=jxx(it)+jyy(it)
!end do
!
!do it=1,ntmax/2+1
!write(301,*)ww(it)/omg1,dlog10(jxx(it))
!write(302,*)ww(it)/omg1,dlog10(jyy(it))
!write(303,*)ww(it)/omg1,dlog10(jtot(it))
!end do
!
!end if
!
!call MPI_FINALIZE(ierr)

end
