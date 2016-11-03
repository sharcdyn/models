program run_hess
!!! Program to calculate the hessian in the models
use molecule_
implicit none
integer nat,npot,ipot

real*8, allocatable :: x(:,:),vel(:,:)
character*2, allocatable :: sat(:)

complex*16, allocatable :: V(:,:),dipole(:,:,:)
logical, allocatable :: egradV(:,:)
complex*16, allocatable :: graddipole(:,:,:,:,:)
complex*16, allocatable :: gradVp(:,:,:,:),gradVm(:,:,:,:)


real*8 :: dh=0.01
real*8 kk
real*8, allocatable :: xmove(:,:),hess(:,:),hess2(:,:)

integer i,j
integer ii,ij,ji,jj

open(1,file="geom.xyz",status="OLD",iostat=i)
if (i.eq.0) then
 write(6,*) "Reading from geom.xyz"
 read(1,*) nat
else
 write(6,*) "Number of atoms"
 read(5,*) nat
endif
allocate(sat(nat),x(nat,3),vel(nat,3))
allocate(xmove(nat,3),hess(3*nat,3*nat),hess2(3*nat,3*nat))

vel=0.
if (i.eq.0) then
 read(1,*) 
 do ii=1,nat
  read(1,*) sat(ii),(x(ii,ij),ij=1,3)
  do ij=1,3
   x(ii,ij)=x(ii,ij)/.5292
  enddo
 enddo
 close(1)
else
 open(2,file="geom")
 do i=1,nat
  read(2,*) sat(i),kk,(x(i,j),j=1,3)
 enddo
 close(2)
endif

write(6,*) "Geometry ready"

write(6,*) "Number of potentials in the model and the potential for the hessian if (negative Hessian will be ignored)"
read(5,*) npot,ipot

allocate (V(npot,npot),egradV(npot,npot),dipole(3,npot,npot))
allocate (gradVp(npot,npot,nat,3),gradVm(npot,npot,nat,3))
allocate (graddipole(3,npot,npot,nat,3))

call modeldescription(nat,npot)

!!! Checking gradient
xmove=x
  call molecule (nat,npot,xmove,sat,V, dipole,gradVp,egradV, graddipole, dh, vel)

write(6,*) "Gradient for the geometry in file grad.dat"
open(1,file="grad.dat")
do i=1,npot
 do j=1,npot
  write(1,*) i,j
  do ii=1,nat
   write(1,900) (dreal(gradVp(i,j,ii,ij)),ij=1,3)
  enddo
 enddo
enddo
close(1)

write(6,*) "Potential matrix in H0.dat"
open(1,file="H0.dat")
do i=1,npot
 write(1,900) (dreal(V(i,j)),j=1,npot)
enddo

write(6,*) "Dipoles in dipole1.dat dipole2.dat dipole3.dat"
do ii=1,3
 if (ii.eq.1) open(1,file="dipole1.dat")
 if (ii.eq.2) open(1,file="dipole2.dat")
 if (ii.eq.3) open(1,file="dipole3.dat")
 do i=1,npot
  write(1,900) (dreal(dipole(ii,i,j)),j=1,npot)
 enddo
enddo
!! Hessian if required

if (ipot.gt.0) then
 i=0
 do ii=1,nat
  do ij=1,3
   i=i+1
   xmove(ii,ij)=xmove(ii,ij)+dh
   call molecule (nat,npot,xmove,sat,V, dipole,gradVp,egradV, graddipole, dh, vel)
   xmove(ii,ij)=xmove(ii,ij)-2*dh
   call molecule (nat,npot,xmove,sat,V, dipole,gradVm,egradV, graddipole, dh, vel)
   xmove(ii,ij)=xmove(ii,ij)+dh
   j=0
   do ji=1,nat
    do jj=1,3
     j=j+1
     hess2(i,j)=dreal(gradVp(ipot,ipot,ji,jj)-gradVm(ipot,ipot,ji,jj))/2./dh
    enddo
   enddo
  enddo
 enddo

 kk=0.
 do i=1,3*nat
  do j=i,3*nat
   kk=kk+(hess2(i,j)-hess2(j,i))**2
   hess(i,j)=(hess2(i,j)+hess2(j,i))/2.
   hess(j,i)=hess(i,j)
  enddo
 enddo
 write(6,*) "Non-Hermitian part of the hessian ",kk
 write(6,*) "Hessian in hessian.dat"
 open(1,file="hessian.dat")
 do i=1,3*nat
  write(1,900) (hess(i,j),j=1,i)
 enddo
endif

900 format(100000(x,E20.10e3))
end program
