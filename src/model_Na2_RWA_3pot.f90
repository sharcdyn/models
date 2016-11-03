module molecule_
implicit none
contains

subroutine modeldescription(nat,npot)
implicit none
integer,intent(in) :: nat,npot

write(6,*) "Model of Na2, taken from Bajo et al. JPCA 116, 2800 (2012)"
write(6,*) "Rotating Wave Approximation with 3 potentials"
write(6,*) "Pot 1 -> V1"
write(6,*) "Pot 2 -> V2-0.08055"
write(6,*) "Pot 3 -> V3-0.08055-0.03606"

if (nat.ne.2) stop 'Potential model for Na2, only two atoms'
if (npot.ne.3) stop 'Model for 3 potentials using RWA'

end

subroutine molecule (nat,npot,x0,sat,V, dipole,gradV,egradV, graddipole, dh, vel)
implicit none
integer, intent(in) :: nat,npot
character*2, intent(in) :: sat(nat)
real*8, intent(in) :: x0(nat,3),vel(nat,3),dh
complex*16, intent(out) :: V(npot,npot),dipole(3,npot,npot)
complex*16,intent(out) :: gradV(npot,npot,nat,3),graddipole(3,npot,npot,nat,3)
logical, intent(out) :: egradV(npot,npot)

real*8 x(nat,3)
complex*16 :: Vp(npot,npot),dipolep(3,npot,npot)
complex*16 :: Vm(npot,npot),dipolem(3,npot,npot)

integer i,j,k
integer ii,ij

!

call potxyz(nat,npot,x0,V,dipole)
!! gradients
x=x0
do ii=1,nat
 do ij=1,3
  x(ii,ij)=x(ii,ij)+dh
  call potxyz(nat,npot,x,Vp,dipolep)
  x(ii,ij)=x(ii,ij)-2.*dh
  call potxyz(nat,npot,x,Vm,dipolem)
  x(ii,ij)=x(ii,ij)+dh
  do i=1,npot
   do j=1,npot
    gradV(i,j,ii,ij)=(Vp(i,j)-Vm(i,j))/2./dh
    egradV(i,j)=.true.
    do k=1,3
     graddipole(k,i,j,ii,ij)=(dipolep(k,i,j)-dipolem(k,i,j))/2./dh
    enddo
   enddo
  enddo
 enddo
enddo
    
end subroutine molecule


subroutine potxyz(nat,npot,x,V,dipole)
implicit none
integer, intent(in) :: nat,npot
real*8, intent(in) :: x(nat,3)
complex*16, intent(out) :: V(npot,npot),dipole(3,npot,npot)

real*8 r

r=sqrt( (x(1,1)-x(2,1))**2+ (x(1,2)-x(2,2))**2+ (x(1,3)-x(2,3))**2 )
!write(6,*) "Distance ",r
call potential(npot,r,V,dipole)

return
end subroutine potxyz

subroutine potential(npot,r,V,dipole)
implicit none
integer, intent(in) :: npot
real*8, intent(in) :: r
complex*16, intent(out) :: dipole(3,npot,npot),V(npot,npot)
!cosas propias
real*8 V0(npot), De(npot), a(npot)
real*8 k(npot),r0(npot),kk
integer i,j,ii

if (npot.ne.3) stop "Only available for 3 potentials"


!pot 1
k(1)= 0.00954206
De(1)= 0.0270664
r0(1)= 5.8249978
V0(1)= 0.0000014
!pot 2
k(2)= 0.00587197
De(2)= 0.0357438
r0(2)= 6.8691376
V0(2)= 0.0665470-0.08055

!pot 3
k(3)= 0.00533354
De(3)= 0.0241634
r0(3)= 6.7425752
V0(3)= 0.1167084-0.08055-0.03606



do ii=1,3
 do i=1,npot
  do j=1,npot
   dipole(ii,i,j)=dcmplx(0.d0,0.d0)
  enddo
 enddo
enddo



kk=3.018 -0.705572*r +0.308706*r**2 -0.039349*r**3+ &
0.0020223*r**4 -3.6789e-05*r**5
if (r.gt.19.99) kk=3.568900063
dipole(1,1,2)=dcmplx(kk,0.d0)

kk=100.276 -99.8748*r +40.0344*r**2 -8.53189*r**3 +1.06116*r**4 &
-0.0797605*r**5 +0.00357366*r**6 -8.7917e-05*r**7 +9.1424e-07*r**8
if (r.gt.19.99) kk=-2.745101531
dipole(2,2,3)=dcmplx(kk,0.d0)

 dipole(1,2,1)=dconjg(dipole(1,1,2))
 dipole(2,3,2)=dconjg(dipole(2,2,3))

do i=1,npot
 do j=1,npot
  V(i,j)=dcmplx(0.,0.)
 enddo
! V(i,i)=dcmplx(k(1,i)*(r(1)-r0(1,i))**2+V0(i),0.d0)
 a(i)=dsqrt(k(i)/(2*De(i)))
 V(i,i)=dcmplx(De(i)*(1.-dexp(-a(i)*(r-r0(i))))**2+V0(i),0.d0)
enddo

return
end subroutine

end module
