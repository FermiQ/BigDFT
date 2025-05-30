!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_symanal
!! NAME
!! abi_symanal
!!
!! FUNCTION
!! Find the space group, Bravais lattice, including Shubnikov characteristics
!! from the list of symmetries (including magnetic characteristics), and lattice parameters
!! Warning : the recognition of the space group might not yet work for the
!! Shubnikov group of type IV
!! 
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (XG, RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! chkprim= if 1 then stop if the cell is not primitive
!! msym=default maximal number of symmetries
!! nsym=actual number of symmetries
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! symafm(1:msym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,1:msym)=symmetry operations in real space in terms
!!  of primitive translations
!! tnons(3,1:msym)=nonsymmorphic translations for symmetry operations
!! tolsym=tolerance for the symmetry operations
!!
!! OUTPUT
!! bravais(11)=characteristics of Bravais lattice (see abi_symlatt.F90)
!! genafm(3)=magnetic translation generator (in case of Shubnikov group type IV)
!! ptgroupma = magnetic point group number
!! spgroup=symmetry space group
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine abi_symanal(bravais,chkprim,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tolsym)

 use abi_defs_basis
 use abi_interfaces_lowlevel
 use abi_interfaces_geometry, except_this_one => abi_symanal

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chkprim,msym,nsym
 integer,intent(out) :: ptgroupma,spgroup
 real(dp),intent(in) :: tolsym
!arrays
 integer,intent(out) :: bravais(11)
 integer,intent(in) :: symafm(msym),symrel(3,3,msym)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: tnons(3,msym)
 real(dp),intent(out) :: genafm(3)

!Local variables-------------------------------
!scalars
 integer :: iholohedry_nomagn,isym,isym_nomagn,multi
 integer :: nptsym,nsym_nomagn,shubnikov
 character(len=5) :: ptgroup,ptgroupha
 character(len=500) :: message
!arrays
 integer :: identity(3,3)
 integer,allocatable :: ptsymrel(:,:,:),symrel_nomagn(:,:,:)
 real(dp),allocatable :: tnons_nomagn(:,:)

! *************************************************************************

!This routine finds the Bravais characteristics, without actually
!looking at the symmetry operations.
 allocate(ptsymrel(3,3,msym))
 call abi_symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
 deallocate(ptsymrel)

!Check whether the cell is primitive or not.
 call abi_chkprimit(chkprim,multi,nsym,symafm,symrel)

 spgroup=0 ; ptgroupma=0 ; genafm(:)=zero
 if(multi>1)then
!  Modify bravais if the cell is not primitive ; no determination of the space group
   bravais(1)=-bravais(1)
 else

!  The cell is primitive, so that the space group can be
!  determined. Need to distinguish Fedorov and Shubnikov groups.
!  Do not distinguish Shubnikov types I and II.
!  Also identify genafm, in case of Shubnikov type IV
   identity(:,:)=reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
   shubnikov=1
   do isym=1,nsym
     if(symafm(isym)==-1)then
       shubnikov=3
       if(sum(abs(symrel(:,:,isym)-identity(:,:)))==0)then
         shubnikov=4
         genafm(:)=tnons(:,isym)
         exit
       end if
     end if
   end do

   if(shubnikov/=1)then
     if(shubnikov==3)write(message, '(a)' )' Shubnikov space group type III'
     if(shubnikov==4)write(message, '(a)' )' Shubnikov space group type IV'
     call abi_wrtout(std_out,message,'COLL')
   end if

   if(shubnikov==1 .or. shubnikov==3)then
!    Find the correct Bravais characteristics and point group
!    Should also be used for Shubnikov groups of type IV ...
     call abi_symbrav(bravais,msym,nsym,ptgroup,rprimd,symrel,tolsym)
!    Find the space group
     call abi_symspgr(bravais,nsym,spgroup,symrel,tnons,tolsym)
   end if

   if(shubnikov/=1)then

!    Determine nonmagnetic symmetry operations
     nsym_nomagn=nsym/2
     allocate(symrel_nomagn(3,3,nsym_nomagn),tnons_nomagn(3,nsym_nomagn))
     isym_nomagn=0
     do isym=1,nsym
       if(symafm(isym)==1)then
         isym_nomagn=isym_nomagn+1
         symrel_nomagn(:,:,isym_nomagn)=symrel(:,:,isym)
         tnons_nomagn(:,isym_nomagn)=tnons(:,isym)
       end if
     end do

     if(shubnikov==3)then

!      Find the point group of the halved symmetry set
       call abi_symptgroup(iholohedry_nomagn,nsym_nomagn,ptgroupha,symrel_nomagn)

!      Deduce the magnetic point group (ptgroupma) from ptgroup and ptgroupha
       call abi_getptgroupma(ptgroup,ptgroupha,ptgroupma)

     else if(shubnikov==4)then

!      Find the Fedorov space group of the halved symmetry set
       call abi_symspgr(bravais,nsym_nomagn,spgroup,symrel_nomagn,tnons_nomagn,tolsym)

!      The magnetic translation generator genafm has already been determined

     end if

     deallocate(symrel_nomagn,tnons_nomagn)    !  added by MM on Oct.25

   end if ! Shubnikov groups

 end if

end subroutine abi_symanal
!!***
