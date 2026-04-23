module variables
implicit none
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14)
! DEBUG NOTE:
! Use a standard double-precision kind selector.
integer,parameter::n=256
real(dp),parameter::rho=0.6360d0,sigma=1.00d0,rc=2.50d0*sigma
real(dp),parameter::boxlength=(real(n,dp)/rho)**(1.0_dp/3.0_dp)
real(dp),parameter::lxh=boxlength/2.0_dp
real(dp),parameter::dt=0.0010d0
! DEBUG FIX:
! The original file used/commented an expression with dfloat(...).
! For portability and standard Fortran, use real literals like 2.0_dp
! and real(n,dp) when integer-to-real conversion is needed.
!(real(n,dp)/rho)**(1.0d0/3.0d0)
real(dp),dimension(n)::x,y,z,x_p,y_p,z_p,vx,vy,vz
real(dp),dimension(n)::xm,ym,zm
real(dp),dimension(n)::fx,fy,fz,f
real(dp)::pot_energy,p_tot,ke,tot_en,ke_tot
integer::i,j,k
real(dp)::xr,yr,zr,dr,fxij,fyij,fzij,fc,ufc,r2
real(dp)::x1,y1,z1,x2,y2,z2
real(dp)::dx,dy,dz,tr
real(dp)::avx,avy,avz,sumv2,vxbar,vybar,vzbar,fs

contains
! DEBUG FIX:
! The original code used nonstandard srand/rand.
! This helper switches to standard-conforming random_seed/random_number.
subroutine init_random_seed(seed)
  integer, intent(in) :: seed
  integer :: nseed, i
  integer, allocatable :: put(:)
  call random_seed(size=nseed)
  allocate(put(nseed))
  do i = 1, nseed
     put(i) = seed + 37*(i-1)
  end do
  call random_seed(put=put)
  deallocate(put)
end subroutine
end module

!****************************************************************
!                MAIN PROGRAM
!****************************************************************
program main
use variables
implicit none
real(dp)::t,tf
t=0.0d0
tf=5.10d0
print*,boxlength
open(unit=76,file='tthundredthlast1.dat',action='write')

call init_pos
call init_vel
do while (t.lt.tf)
call force_calc
call integrate

tot_en=ke_tot+p_tot
print*,'total energy::',tot_en,'kinetic',ke_tot,'potential::',p_tot
write(76,*)t,tot_en,p_tot,ke_tot
t=t+dt
end do
end program
!********************************************************************
!                INITIALISATON OF POSITION
!********************************************************************
subroutine init_pos
use variables 
implicit none

! DEBUG NOTE:
! Fixed seed for reproducibility while debugging.
call init_random_seed(34)
call random_number(x_p(1))
x_p(1)=x_p(1)*boxlength
call random_number(z_p(1))
z_p(1)=z_p(1)*boxlength
call random_number(y_p(1))
y_p(1)=y_p(1)*boxlength
print*,x_p(1),y_p(1),z_p(1)
do i=2,n
20	call random_number(x_p(i))
	x_p(i)=x_p(i)*boxlength
	call random_number(y_p(i))
	y_p(i)=y_p(i)*boxlength
	call random_number(z_p(i))
	z_p(i)=z_p(i)*boxlength
	do j=1,i-1
		dx=x_p(i)-x_p(j)
		dy=y_p(i)-y_p(j)
		dz=z_p(i)-z_p(j)
		tr=sqrt(dx*dx+dy*dy+dz*dz)
		! WARNING:
		! This overlap check does NOT use minimum-image periodic distance.
		! Two particles near opposite box faces may still be close under PBC.
		
		if (tr.lt. rc) then
			goto 20
		
		end if
	print*,tr,'boxlength',boxlength
	end do
end do
do k=1,n
print*,k,'x::',x_p(k),'y::',y_p(k),'z::',z_p(k)
end do

end subroutine

!***************************************************************
!             INITIALISATION OF VELOCITIES
!***************************************************************
subroutine init_vel

use variables
implicit none


! DEBUG FIX:
! Replaced legacy rand() calls with standard random_number().
do i=1,n
	call random_number(vx(i))
	vx(i)=(vx(i)-0.50d0)*2.075d0
	call random_number(vy(i))
	vy(i)=(vy(i)-0.50d0)*2.075d0
	call random_number(vz(i))
	vz(i)=(vz(i)-0.50d0)*2.075d0
end do

do i=1,n
	x(i)=x_p(i)+vx(i)*dt
	y(i)=y_p(i)+vy(i)*dt
	z(i)=z_p(i)+vz(i)*dt
	print*,'initial:',x_p(i),'calculated:',x(i)
end do

end subroutine
!********************************************************************
!             FORCE CALCULATION
!********************************************************************
subroutine force_calc
use variables
implicit none
p_tot=0.0d0
fc=48.0d0*rc*((sigma/rc)**14-0.50d0*(sigma/rc)**8)
ufc=fc*rc+4.0d0*((sigma/rc)**12-(sigma/rc)**6)
! WARNING:
! fc and ufc are computed, but they are not actually applied below.
! So the code is NOT using a shifted/force-corrected cutoff yet.

do i=1,n 	!setting forces to zero
	fx(i)=0
	fy(i)=0
	fz(i)=0
end do
!let's begin force and potential calculation
do i=1,n-1
	do j=i+1,n
		xr=x(i)-x(j)
		yr=y(i)-y(j)
		zr=z(i)-z(j)
		if (abs(xr).ge.lxh) xr=(boxlength-abs(xr))*((-1.0d0*xr)/abs(xr))
		if (abs(yr).ge.lxh) yr=(boxlength-abs(yr))*((-1.0d0*yr)/abs(yr))
		if (abs(zr).ge.lxh) zr=(boxlength-abs(zr))*((-1.0d0*zr)/abs(zr))
		dr=sqrt(xr*xr+yr*yr+zr*zr)
		! WARNING:
		! The original cutoff logic is effectively disabled because the if-block
		! is commented out. That means all pairs interact, regardless of rc.
		!if (dr.lt.rc) then
			pot_energy=4.0d0*((sigma/dr)**12-(sigma/dr)**6)
			fxij=48.0d0*xr*((sigma/dr)**14-0.50d0*(sigma/dr)**8)
			fyij=48.0d0*yr*((sigma/dr)**14-0.50d0*(sigma/dr)**8)
			fzij=48.0d0*zr*((sigma/dr)**14-0.50d0*(sigma/dr)**8)
			fx(i)=fx(i)+fxij
			fy(i)=fy(i)+fyij
			fz(i)=fz(i)+fzij
		
			fx(j)=fx(j)-fxij
			fy(j)=fy(j)-fyij
			fz(j)=fz(j)-fzij
			p_tot=p_tot+pot_energy
		!end if
	end do !j loop
end do !i loop
print*,'potential energy',p_tot
end subroutine
!*****************************************************************************
!                     INTEGRATION OF MOTION
!*****************************************************************************
SUBROUTINE integrate
use variables
implicit none
ke=0.0d0
ke_tot=0.0d0
print*,'refresshing energies',ke,ke_tot

! WARNING:
! These resets are unnecessary because vx, vy, vz are overwritten below.
! They do not break compilation, but removing them makes debugging clearer.
do i=1,n !setting velocities to zero
	vx(i)=0
	vy(i)=0
	vz(i)=0
end do

do i=1,n
	! Verlet-style position update
	xm(i)=2.0d0*x(i)-x_p(i)+fx(i)*(dt*dt)
	ym(i)=2.0d0*y(i)-y_p(i)+fy(i)*(dt*dt)
	zm(i)=2.0d0*z(i)-z_p(i)+fz(i)*(dt*dt)
	!print*,'integrationof positions',xm(i),ym(i),zm(i)
	vx(i)=(xm(i)-x_p(i))/(2.0d0*dt)
	vy(i)=(ym(i)-y_p(i))/(2.0d0*dt)
	vz(i)=(zm(i)-z_p(i))/(2.0d0*dt)
	!print*,'caluclated velocitites',vx(i),vy(i),vz(i)
	ke_tot=ke_tot+0.50d0*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
	!update old coordinates
	x_p(i)=x(i)
	y_p(i)=y(i)
	z_p(i)=z(i)
	!update new coordinates
	x(i)=xm(i)
	y(i)=ym(i)
	z(i)=zm(i)
	!putting the particles back into box
	if (x(i).gt.boxlength) then
		x(i)=x(i)-boxlength
		x_p(i)=x_p(i)-boxlength
	elseif (x(i).lt.(0.0d0)) then
		x(i)=x(i)+boxlength
		x_p(i)=x_p(i)+boxlength
	endif
	if (y(i).gt.boxlength) then
		y(i)=y(i)-boxlength
		y_p(i)=y_p(i)-boxlength
	elseif (y(i).lt.(0.0d0)) then
		y(i)=y(i)+boxlength
		y_p(i)=y_p(i)+boxlength
	endif
	if (z(i).gt.boxlength) then
		z(i)=z(i)-boxlength
		z_p(i)=z_p(i)-boxlength
	elseif (z(i).lt.(0.0d0)) then
		z(i)=z(i)+boxlength
		z_p(i)=z_p(i)+boxlength
	endif

end do
print*,'kinetic energy total',ke_tot
end subroutine

