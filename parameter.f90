module parameter
implicit none
include'mpif.h'
include'fftw3.f'
integer(kind=4),parameter::nvar=6
real(8),parameter::pi=3.141592653589793D+00,ppi=pi
real(8),parameter::vf=(1*1.d6)/(2.1876879*1.d6)
!real(8),parameter::omg1=45.56d0/3200d0
!real(8),parameter::period=2*pi/omg1
!real(8),parameter::esi=0,tao=3.0*period,cep=0.5*ppi
!real(8),parameter::fullperiod=6*tao*(1+esi*esi)**0.5d0
!integer,parameter::ntmax=1024*4*16,ntmsax=1024*8*16
!real(8),parameter::dt=fullperiod/dble(ntmax),dts=fullperiod/dble(ntmsax)
double complex,parameter::ur=(1.d0,0.d0),ui=(0.d0,1.d0)
integer,parameter::ntm=1024*16,ntmax=ntm*2,ntms=1024*32,ntmsax=ntms*2
real(8),parameter::fm1=(2.75*1.d6)/(5.1422062*1.d11),omg1=1.d0*2*ppi/(4.134138*1.d4)
real(8),parameter::aa0=fm1/omg1,period=2*ppi/omg1,pul=4.2
real(8),parameter::dt=period/dble(ntm),dts=period/dble(ntms)
!real(8),parameter::I0=5.0*1.d11,iau=3.51*1.d16
!real(8),parameter::fm1=(I0/Iau/(1+esi*esi)**0.5d0)**0.5d0
!real(8),parameter::aa0=fm1/omg1
integer,parameter::nnk=40,nnt=40
real(8),parameter::tss1=2.5d0/period,epa=1.d0

end module
