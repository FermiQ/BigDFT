!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_xchelu
!! NAME
!! abi_xchelu
!!
!! FUNCTION
!! Returns exc, vxc, and eventually d(vxc)/d($\rho$) from input rho.
!!
!! NOTES
!! Hedin-Lundqvist exchange and correlation (xc)--
!! L. Hedin and B.I. Lundqvist, J. Phys. C. 4, 2064 (1971).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR, LG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npt=number of real space points on which density is provided
!!  order=gives the maximal derivative of Exc computed.
!!  rspts(npt)=Wigner-Seitz radii at each point
!!
!! OUTPUT
!!  exc(npt)=exchange-correlation energy density (hartree)
!!  vxc(npt)=xc potential (d($\rho$*exc)/d($\rho$)) (hartree)
!!  if(order>1) dvxc(npt)=derivative d(vxc)/d($\rho$) (hartree*bohr^3)
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine abi_xchelu(exc,npt,order,rspts,vxc,dvxc)  ! dvxc is optional

 use abi_defs_basis
 use abi_interfaces_lowlevel

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npt,order
!arrays
 real(dp),intent(in) :: rspts(npt)
 real(dp),intent(out) :: exc(npt),vxc(npt)
 real(dp),intent(out),optional :: dvxc(npt)

!Local variables-------------------------------
!aa and cc are H-L fitting parameters A and C (C in hartree)
!rs = (3/(4 Pi))**(1/3) * rho(r)**(-1/3).
!scalars
 integer :: ipt
 real(dp),parameter :: aa=21_dp,c1_21=one/21_dp,c4_9=4.0_dp/9.0_dp,cc=0.0225_dp
 real(dp) :: dfac,efac,rs,rsm1,vfac,xx
 character(len=500) :: message

! *************************************************************************

!Checks the values of order
 if(order<0 .or. order>2)then
   write(message, '(a,a,a,a,a,a,i6,a)' )ch10,&
&   ' abi_xchelu : BUG -',ch10,&
&   '  With Hedin-Lundqvist xc functional, the only',ch10,&
&   '  allowed values for order are 0, 1 or 2, while it is found to be',&
&   order,'.'
   call abi_wrtout(std_out,message,'COLL')
   call abi_leave_new('COLL')
 end if

!Compute vfac=(3/(2*Pi))^(2/3)
 vfac=(1.5_dp/pi)**(2.0_dp/3.0_dp)
!Compute efac=(3/4)*vfac
 efac=0.75_dp*vfac
!Compute dfac=(4*Pi/9)*vfac
 dfac=(4.0_dp*pi/9.0_dp)*vfac
!separate cases with respect to order
 if (order==2) then
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=one/rs
!    compute energy density exc (hartree)
     xx=rs*c1_21
     exc(ipt)=-cc*((one+xx**3)*log(one+one/xx)+&
&     half*xx-xx*xx-third) - efac*rsm1
!    compute xc potential d(rho*exc)/d(rho) (hartree)
     vxc(ipt)=-cc*log(one+aa*rsm1)-vfac*rsm1
!    compute d(vxc)/d(rho) (hartree*bohr^3)
     dvxc(ipt)=-(rs**2)*((c4_9*pi)*cc*rs/(one+xx) + dfac)
   end do
 else
!  Loop over grid points
   do ipt=1,npt
     rs=rspts(ipt)
     rsm1=one/rs
!    compute energy density exc (hartree)
     xx=rs*c1_21
     exc(ipt)=-cc*((one+xx**3)*log(one+one/xx)+&
&     half*xx-xx*xx-third) - efac*rsm1
!    compute xc potential d(rho*exc)/d(rho) (hartree)
     vxc(ipt)=-cc*log(one+aa*rsm1)-vfac*rsm1
   end do
 end if

end subroutine abi_xchelu
!!***
