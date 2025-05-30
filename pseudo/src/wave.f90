!> @file
!! atomic program for generating and optimizing hgh pseudo-potentials.
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by matthias krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/
!!
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Calc. int(psi(r)*dr)
real(kind=8) function Wwav(ng,ll,xp,psi,r)
   implicit none
   !Arguments
   integer, intent(in) :: ng,ll
   real(kind=8), dimension(0:ng), intent(in) :: psi,xp
   real(kind=8), intent(in) :: r
   !Local variables
   integer :: i
   real(kind=8) :: t1,t13,t19,t2,t3,t4,t8,t9
   wwav=0.d0
   if (ll.eq.0) then
      t1=sqrt(0.3141592653589793D1)/2.0d0
      do i=0,ng
         t3 = sqrt(xp(i))
         wwav=wwav+1.D0/t3*t1*Derf(t3*r)*psi(i)
      enddo
   elseif (ll.eq.1) then
      t2 = r*r
      do i=0,ng
         t1 = xp(i)
         t9 = -exp(-t1*t2)/t1*psi(i)*0.5d0
         wwav=wwav+t9
      enddo
   elseif (ll.eq.2) then
      t2 = r*r
      t3=sqrt(0.3141592653589793D1)/4.D0
      do i=0,ng
         t1 = xp(i)
         t8 = sqrt(t1)
         t9 = t8*t8
         t19 = (-r*exp(-t1*t2)/t1/2.D0+1.D0/t9/t8 &
              *t3*Derf(t8*r))*psi(i)
         wwav=wwav+t19
      enddo
   elseif (ll.eq.3) then
      t1 = r*r
      do i=0,ng
         t2 = xp(i)
         t4 = exp(-t1*t2)
         t8 = t2*t2
         t13 = (-t1*t4/t2/2.D0-1.D0/t8*t4/2.D0)*psi(i)
         wwav=wwav+t13
      enddo
   else
      stop 'sorry l>3 not implemented'
   endif
end function Wwav


!> Calculates psi(r)
real(kind=8) function wave(ng,ll,xp,psi,r)
   implicit none
   !Arguments
   integer, intent(in) :: ng,ll
   real(kind=8), dimension(0:ng), intent(in) :: psi,xp
   real(kind=8), intent(in) :: r
   !Local variables
   integer :: i
   real(kind=8) :: t3

   wave=0.d0
   t3=r*r
   if (ll.eq.0) then
      do i=0,ng
         wave=wave + psi(i)*exp(-xp(i)*t3)
      enddo
   elseif (ll.eq.1) then
      do i=0,ng
           wave=wave + r*psi(i)*exp(-xp(i)*t3)
      enddo
   elseif (ll.eq.2) then
      do i=0,ng
         wave=wave + t3*psi(i)*exp(-xp(i)*t3)
      enddo
   elseif (ll.eq.3) then
      do i=0,ng
         wave=wave + t3*r*psi(i)*exp(-t3*xp(i))
      enddo
   else
      stop 'sorry l>3 not implemented'
   endif
   
end function wave


!> Calculates first derivative of psi(r)
real(kind=8) function dwave(ng,ll,xp,psi,r)
   implicit none
   !Arguments
   integer, intent(in) :: ng,ll
   real(kind=8), dimension(0:ng), intent(in) :: psi,xp
   real(kind=8), intent(in) :: r
   !Local variables
   integer :: i
   real(kind=8) :: t1,t2,t3,t4,t5,t6,t8,t9,t10,t12

   dwave=0.d0
   t4 = r*r
   if (ll.eq.0) then
      do i=0,ng
         t2 = xp(i)
         t9 = -2.D0*psi(i)*t2*r*exp(-t2*t4)
         dwave=dwave + t9
      enddo
   elseif (ll.eq.1) then
      do i=0,ng
         t1 = psi(i)
         t2 = xp(i)
         t5 = exp(-t2*t4)
         t10 = t1*t5-2.D0*t4*t1*t2*t5
         dwave=dwave + t10
      enddo
   elseif (ll.eq.2) then
      do i=0,ng
         t1 = psi(i)
         t3 = xp(i)
         t6 = exp(-t3*t4)
         t12 = 2.D0*r*t1*t6-2.D0*t4*r*t1*t3*t6
         dwave=dwave + t12
      enddo
   elseif (ll.eq.3) then
      do i=0,ng
         t2 = psi(i)
         t5 = xp(i)
         t6 = exp(-t5*t4)
         t8 = t4**2
         t12 = 3.D0*t4*t2*t6-2.D0*t8*t2*t5*t6
         dwave=dwave + t12
      enddo
   else
      stop 'sorry l>3 not implemented'
   endif
   
end function dwave


!> Calculates second derivative of psi(r)
real(kind=8) function ddwave(ng,ll,xp,psi,r)
   implicit none
   !Arguments
   integer, intent(in) :: ng,ll
   real(kind=8), dimension(0:ng), intent(in) :: psi,xp
   real(kind=8), intent(in) :: r
   !Local variables
   integer :: i
   real(kind=8) :: t1,t2,t3,t4,t5,t6,t8,t10,t11,t12,t14,t15,t18

   ddwave=0.d0
   t4=r*r
   if (ll.eq.0) then
      do i=0,ng
         t1 = psi(i)
         t2 = xp(i)
         t6 = exp(-t2*t4)
         t8 = t2*t2
         t12 = -2.D0*t1*t2*t6+4.D0*t1*t8*t4*t6
         ddwave=ddwave+t12
      enddo
   elseif (ll.eq.1) then
      do i=0,ng
         t1 = psi(i)
         t2 = xp(i)
         t6 = exp(-t2*t4)
         t11 = t2**2
         t14 = -6.D0*t1*t2*r*t6+4.D0*t4*r*t1*t11*t6
         ddwave=ddwave+t14
      enddo
   elseif (ll.eq.2) then
      do i=0,ng
         t1 = psi(i)
         t2 = xp(i)
         t5 = exp(-t2*t4)
         t10 = t4*t4
         t12 = t2*t2
         t15 = 2.D0*t1*t5-10.D0*t4*t1*t2*t5+4.D0*t10*t1*t12*t5
         ddwave=ddwave+t15
      enddo
   elseif (ll.eq.3) then
      do i=0,ng
         t1 = psi(i)
         t3 = xp(i)
         t6 = exp(-t3*t4)
         t12 = t4*t4
         t15 = t3*t3
         t18 = 6.D0*r*t1*t6-14.D0*t4*r*t1*t3*t6+4.D0*t12*r*t1 &
              *t15*t6
         ddwave=ddwave+t18
      enddo
   else
      stop 'sorry l>3 not implemented'
   endif

end function ddwave
