module molecule_
implicit none
contains

subroutine modeldescription(nat,npot)
implicit none
integer, intent(in) :: nat,npot

write(6,*) "Model of IBr taken from Marquetand et al. Faraday Disc. 153, 261 (2011)"
write(6,*) "Potential part is from Guo JCP 99, 1685 (1993)"
write(6,*) "Dipole matrix adjusted from Patchkovskii PCCP 6, 926 (2006)"

if (nat.ne.2) stop "Number of atoms is not valid, must be 2"
if (npot.ne.3) stop "Numer of potentials is not valid, must be 3"

end subroutine

subroutine molecule (nat,npot,x,sat,V, dipole,gradV,egradV, graddipole, dh, vel)
implicit none
integer, intent(in) :: npot,nat
real*8, intent(in) ::  x(nat,3),dh, vel(nat,3)
character*2, intent(in) :: sat(nat)

complex*16, intent(out) :: V(npot,npot), dipole(3,npot,npot)
complex*16, intent(out) :: graddipole(3,npot,npot,nat,3), gradV(npot,npot,nat,3)
logical, intent(out) :: egradV(npot,npot)

integer i,j


egradV=.true.

call potential( nat,npot,x,V)
call gradpot( nat,npot,x,dh,gradV)

call dip( nat,npot,x, dipole)
call graddip( nat, npot, x, dh, graddipole)

return
end

subroutine potential( nat,npot,x,V )
implicit none
integer, intent(in) :: nat,npot
real*8, intent(in) :: x(nat,3)
complex*16, intent(out) :: V(npot,npot)

real*8 r
integer i

r=0.
do i=1,3
 r=r+(x(1,i)-x(2,i))**2
enddo
r=sqrt(r)

call potr(npot,r,V)

return
end subroutine potential

subroutine potr(npot,r,V)
implicit none
integer, intent(in) :: npot
real*8, intent(in) :: r
complex*16, intent(out) :: V(npot,npot)

real*8 V0,V1,V2

if (npot.ne.3) stop "Numer of potentials is not valid, must be 3"

V0=0.067*((1-exp(-0.996*(r-4.666)))**2-1)
V1=0.01019*((1-exp(-1.271*(r-5.3479)))**2-1)+0.01679
V2=2.82*exp(-0.9186*r)+3.0e7*exp(-4.3*r)

V=cmplx(0.,0.)
V(1,1)=V0*(1.,0.)
V(2,2)=V1*(1.,0.)
V(3,3)=V2*(1.,0.)
V(2,3)=0.6834e-3*(1.,0.)
V(3,2)=dconjg(V(2,3))

return
end subroutine potr

subroutine gradpot(nat,npot,x,dh,gradV)
implicit none
integer, intent(in) :: nat,npot
real*8, intent(in) :: x(nat,3)
complex*16, intent(out) :: gradV(npot,npot,nat,3)
real*8, intent(in) :: dh

real*8 xmove(nat,3)
complex*16 Vm(npot,npot),Vp(npot,npot)
integer i,j

xmove=x
do i=1,nat
 do j=1,3
  xmove(i,j)=xmove(i,j)+dh
  call potential( nat,npot,xmove,Vp)
  xmove(i,j)=xmove(i,j)-2*dh
  call potential( nat,npot,xmove,Vm)
  xmove(i,j)=xmove(i,j)+dh
  gradV(:,:,i,j)=(Vp-Vm)/2./dh
 enddo
enddo

return
end subroutine gradpot


subroutine dip( nat,npot,x, dipole)
implicit none
integer, intent(in) :: nat,npot
real*8, intent(in) :: x(nat,3)
complex*16, intent(out) :: dipole(3,npot,npot)

integer i
real*8 r

if (nat.ne.2) stop "Number of atoms is not valid, must be 2"
if (npot.ne.3) stop "Numer of potentials is not valid, must be 3"

r=0.
do i=1,3
 r=r+(x(1,i)-x(2,i))**2
enddo
r=sqrt(r)

call dipr(npot,r,dipole)

return
end subroutine dip

subroutine dipr( npot,r, dipole)
implicit none
integer, intent(in) :: npot
real*8, intent(in) :: r
complex*16, intent(out) :: dipole(3,npot,npot)

real*8 :: A=1./.5292,D=1./2.541746230211

real*8 a4, a10 , a11 , a12 , a13 , a14 , a15 , a16 , a17 , a18 , a19 , a20 
real*8 a21 , a22 , a23 , a24 , a25 , a26 , a27 , a28 , a29 , a30 , a31 
real*8 a32 , a33 , a34 , a35 , a36 , a37
a4= 8.20 *A**3
a10=2.95 *A
a11=50.02*A**3
a12=5.62 /A**2
a13=1.68 *A
a14=3.15 *A**3
a15=3.5  /A**2
a16=3.13 *A
a17=2.25 *D/A
a18=2.07 *A
a19=0.83 /A**2
a20=0.80 *D
a21=4.00 /A
a22=2.25 *A
a23=0.30 *D
a24=20.00/A**2
a25=2.85 *A
a26=0.90 *D
a27=5.00 /A**2
a28=2.65 *A
a29=0.45 *D
a30=5.00 /A**2
a31=3.35 *A
a32=0.59 *D/A
a33=2.07 *A
a34=0.67 /A**2
a35=0.38 *D/A
a36=2.41 *A
a37=3.58 /A**2

if (npot.ne.3) stop "Numer of potentials is not valid, must be 3"


dipole=(0.,0.)

dipole(3,1,1)=a17*(r-a18)*exp(-a19*(r-a18)**2)
dipole(3,2,2)=a20*(1-exp(-a21*(r-a22)**2))-a20+a23*exp(-a24*(r-a25)**2)
dipole(3,3,3)=a26*exp(-a27*(r-a28)**2)+a29*exp(-a30*(r-a31)**2)
dipole(3,1,2)=a32*(R-a33)*exp(-a34*(r-a33))
dipole(3,1,3)=a35*(r-a36)*exp(-a37*(r-a36)**2)
dipole(3,3,1)=conjg(dipole(3,1,3))
dipole(3,2,1)=conjg(dipole(3,1,2))

return
end subroutine dipr

subroutine graddip( nat, npot, x, dh, graddipole)
implicit none
integer, intent(in) :: nat,npot
real*8, intent(in) :: x(nat,3)
complex*16, intent(out) :: graddipole(3,npot,npot,nat,3)
real*8, intent(in) :: dh

real*8 xmove(nat,3)
complex*16 dipm(3,npot,npot),dipp(3,npot,npot)
integer i,j,k

xmove=x
do i=1,nat
 do j=1,3
  xmove(i,j)=xmove(i,j)+dh
  call dip( nat,npot,xmove,dipp)
  xmove(i,j)=xmove(i,j)-2*dh
  call dip( nat,npot,xmove,dipm)
  xmove(i,j)=xmove(i,j)+dh
  graddipole(:,:,:,i,j)=(dipp-dipm)/2./dh
 enddo
enddo

return
end subroutine graddip

end module
