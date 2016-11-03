module molecule_
implicit none

PUBLIC :: paramnm

type T_paramnm
  integer nat, npot
  integer, allocatable :: nm
  real*8, allocatable :: x0(:,:,:), dx(:,:,:,:),k(:,:),E0(:)
  
  contains
    procedure :: init

end type T_paramnm

type(T_paramnm) :: paramnm

contains

subroutine init(this,nat,npot)
implicit none
class(T_paramnm), intent(out) :: this
integer, intent(in) :: nat,npot
integer i,j,k,ipot


character*100 at
character*100 fich
real*8 kk,kk2
real*8, allocatable :: mass(:),massred(:)

this%npot=npot
this%nat=nat

allocate(this%E0(npot))
do ipot=1,npot
 if (npot.eq.1) then
  this%E0(npot)=0.
  write(fich,"(A)") "freqs.xyz"
 else
  write(at,"(I10)") ipot
  write(fich,"(A,A,A)") "freqs",trim(at),".xyz"
 endif
 write(6,*) "Reading normal modes from file ",trim(fich)
 open(1,file=fich,status="old",iostat=i)
 if (i.ne.0) then
  write(6,*) "There is no file for frequencies ",trim(fich)
  stop
 endif
 if (ipot.eq.1) then
  this%nm=0
  write(6,*) "Reading the number of normal modes, must be the same in all potentials"
  do
!!!! if there is more than one potential the energy E0 must be here
   if (npot.ne.1) read(1,*) kk
   read(1,*,end=1) i
   this%nm=this%nm+1
   if (i.ne.nat) then
    write(6,*) "Number of atoms in the input ",nat
    write(6,*) "Number of atoms in the frequency calculation ",i
    stop
   endif
   read(1,*)
   do i=1,this%nat
    read(1,*) at,(kk,j=1,3)
   enddo
  enddo

1 continue
  write(6,*) "Number of normal modes ",this%nm
  rewind(1)
  allocate( this%x0(npot,this%nat,3), this%k(npot,this%nm) )
  allocate( this%dx(npot,this%nm,this%nat,3) )
 endif

 do i=1,this%nm
  if (npot.ne.1) read(1,*) this%E0(ipot)
  read(1,*) j
  read(1,*) this%k(ipot,i)
  do j=1,this%nat
   read(1,*) at,(this%x0(ipot,j,k),k=1,3),(this%dx(ipot,i,j,k),k=1,3)
  enddo
 enddo
 close(1)

 do i=1,this%nat
  do j=1,3
   this%x0(ipot,i,j)=this%x0(ipot,i,j)/.5292
   do k=1,this%nm
    this%dx(ipot,k,i,j)=this%dx(ipot,k,i,j)/.5292
   enddo
  enddo
 enddo

 if (ipot.eq.1) then
  allocate(mass(this%nat),massred(this%nm))
  write(6,*) "Reading the masses from geom"
  open(1,file="geom")
  do i=1,this%nat
   read(1,*) at,kk,(kk,j=1,3),mass(i)
   mass(i)=mass(i)*1822.888
  enddo
  close(1)
 endif

 write(6,*) "Potential ",ipot
 do i=1,this%nm
  kk=this%k(ipot,i)
  massred(i)=0.
  do j=1,this%nat
   do k=1,3
    massred(i)=massred(i)+( this%dx(ipot,i,j,k)**2/mass(j))
   enddo
  enddo
  massred(i)=1./massred(i)
  if (this%k(ipot,i).ge.0) then
   this%k(ipot,i)=+ kk*kk*massred(i)
  else
   write(6,*) "Negative frecuency in normal mode ",i
   write(6,*) "Considering a saddle point"
   this%k(ipot,i)=- kk*kk*massred(i)
  endif
  write(6,*) "Normal mode ",i
  write(6,*) "Frequency ",kk,kk*2.194746E5
  write(6,*) "Mass ",massred(i),massred(i)/1822.888
  write(6,*) "Force constant ",this%k(ipot,i)
 enddo
 

 do i=1,this%nm
  kk=0.
  do j=1,this%nat
   do k=1,3
    this%dx(ipot,i,j,k)=this%dx(ipot,i,j,k)*1
    kk=kk+this%dx(ipot,i,j,k)**2
   enddo
  enddo
  kk=sqrt(kk)
!  write(6,*) "Renormalizing displacement for normal mode ",i,kk
  do j=1,this%nat
   do k=1,3
    this%dx(ipot,i,j,k)=this%dx(ipot,i,j,k)/kk*sqrt(mass(j)/massred(i))
   enddo
  enddo
 enddo

enddo

end subroutine

subroutine modeldescription(nat,npot)
implicit none
integer, intent(in) :: nat,npot

write(6,*) "Harmonic potential from generated initial conditions"
call paramnm%init(nat,npot)
write(6,*) "Parameters ready"

end subroutine


subroutine molecule (nat,npot,x,sat,V, dipole,gradV, egradV, graddipole, dh, vel)
implicit none
integer, intent(in) :: npot,nat
real*8, intent(in) :: x(nat,3),vel(nat,3),dh
character*2, intent(in) :: sat(nat)

complex*16, intent(out) :: V(npot,npot), dipole(3,npot,npot)
complex*16, intent(out) :: graddipole(3,npot,npot,nat,3), gradV(npot,npot,nat,3)
logical, intent(out) :: egradV(npot,npot)

real*8 kk, diffnm
real*8 xcopy(nat,3)

integer i,j,k,l,m

if (.not.allocated(paramnm%x0)) then
  write(*,*) "Paramnm was not initialized correctly!"
  stop
endif


do l=1,npot
 V(l,l)=energynm(nat,x,l)
 egradV(l,l)=.true.
enddo
xcopy=x
do i=1,nat
 do j=1,3
  xcopy(i,j)=xcopy(i,j)+dh
  do l=1,npot
   gradV(l,l,i,j)=energynm(nat,xcopy,l)
  enddo
  xcopy(i,j)=xcopy(i,j)-2*dh
  do l=1,npot
   gradV(l,l,i,j)=(gradV(l,l,i,j)-energynm(nat,xcopy,l))/2./dh
  enddo
  xcopy(i,j)=xcopy(i,j)+dh
 enddo
enddo

do i=1,3
 do j=1,npot
  do k=1,npot
   dipole(i,j,k)=dcmplx(0.d0,0.d0)
  enddo
 enddo
enddo

do i=1,npot
 do j=1,npot
  do k=1,nat
   do l=1,3
    do m=1,3
     graddipole(m,i,j,k,l)=dcmplx(0.d0,0.d0)
    enddo
   enddo
  enddo
 enddo
enddo

return
end

real*8 function energynm(nat,x,ipot)
  implicit none
  integer, intent(in) :: nat,ipot
  real*8, intent(in) :: x(nat,3)
  integer i, j, k
  real*8 diffnm

  energynm=paramnm%E0(ipot)
  do i=1, paramnm%nm
    diffnm=0.
    do j=1,nat
      do k=1, 3
        diffnm=diffnm+(x(j,k)-paramnm%x0(ipot,j,k))*paramnm%dx(ipot,i,j,k)
      enddo
    enddo
!    write(6,*) "Movement NM ",i,sqrt(diffnm),paramnm%k(i)*diffnm/2
    energynm=energynm+paramnm%k(ipot,i)*diffnm**2/2
  enddo
  return
end function energynm

end module
