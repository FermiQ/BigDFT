!{\src2tex{textfont=tt}}
!!****f* ABINIT/abi_pawmknhat
!! NAME
!! abi_pawmknhat
!!
!! FUNCTION
!! PAW only:
!! Compute compensation charge density (and derivatives) on the fine FFT grid
!! Can also compute first-order compensation charge density (RF calculations)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  distribfft<type(distribfft_type)>=--optional-- contains infos related to FFT parallelism
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!  ider= 0: nhat(r) is computed
!!        1: cartesian derivatives of nhat(r) are computed
!!        2: nhat(r) and derivatives are computed
!!  idir=direction of atomic displacement (in case of atomic displ. perturb.)
!!  ipert=index of perturbation; must be 0 for ground-state calculations
!!  izero=if 1, unbalanced components of nhat(g) have to be set to zero
!!  me_g0=--optional-- 1 if the current process treat the g=0 plane-wave (only needed when mpi_comm_fft is present)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  mpi_comm_atom=--optional-- MPI communicator over atoms
!!  mpi_comm_fft=--optional-- MPI communicator over FFT components
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  nfft=number of point on the rectangular fft grid
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhatgrdim= -PAW only- 0 if pawgrnhat array is not used ; 1 otherwise
!!  ntypat=number of types of atoms in unit cell.
!!  paral_kgb=--optional-- 1 if "band-FFT" parallelism is activated (only needed when mpi_comm_fft is present)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(my_natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                         (1st-order occupancies if ipert>0)
!!  pawrhoij0(my_natom) <type(pawrhoij_type)>= GS paw rhoij occupancies and related data (used only if ipert>0)
!!                                          set equat to pawrhoij for GS calculations
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  qphon(3)=wavevector of the phonon (RF only)
!!  rprimd(3,3)=dimensional primitive translations for real space
!!  ucvol=volume of the unit cell
!!  xred(3,natom)= reduced atomic coordinates
!!
!! OUTPUT
!!  === if ider=0 or 2
!!    compch_fft=compensation charge inside spheres computed over fine fft grid
!!    pawnhat(nfft,ispden)=nhat on fine rectangular grid
!!  === if ider=1 or 2
!!    pawgrnhat(nfft,ispden,3)=derivatives of nhat on fine rectangular grid (and derivatives)
!!
!! PARENTS
!!      bethe_salpeter,energy,nres2vres,odamix,paw_qpscgw,abi_pawmkrho,respfn,scfcv
!!      scfcv3,screening,setup_positron,sigma
!!
!! CHILDREN
!!      destroy_distribfft,fourdp,free_my_atmtab,get_my_atmtab
!!      init_distribfft_seq,initmpi_seq,abi_mean_fftr,pawexpiqr,pawgylm,pawnhatfr
!!      set_mpi_enreg_fft,timab,unset_mpi_enreg_fft,xmpi_sum,zerosym
!!
!! SOURCE

#include "../libpaw/libpaw.h"

subroutine abi_pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,&
&          my_natom,natom,nfft,ngfft,nhatgrdim,nspden,ntypat,pawang,pawfgrtab,&
&          pawgrnhat,pawnhat,pawrhoij,pawrhoij0,pawtab,qphon,rprimd,ucvol,usewvl,xred,&
&          mpi_atmtab,mpi_comm_atom,mpi_comm_fft,mpi_comm_wvl,me_g0,paral_kgb,distribfft) ! optional arguments

 USE_DEFS
 USE_MSG_HANDLING
 USE_MEMORY_PROFILING
 USE_MPI_WRAPPERS

 !use defs_abitypes,  only : mpi_type
 !use m_mpinfo,       only : set_mpi_enreg_fft, unset_mpi_enreg_fft

 use m_abi_distribfft,   only : distribfft_type, init_distribfft_seq, destroy_distribfft

 use m_pawang,       only : pawang_type
 use m_pawtab,       only : pawtab_type
 use m_pawfgrtab,    only : pawfgrtab_type
 use m_pawrhoij,     only : pawrhoij_type
 use m_paw_finegrid, only : pawgylm, pawexpiqr
 use m_paral_atom,   only : get_my_atmtab, free_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_pawmknhat'
use abi_interfaces_add_libpaw, only: abi_pawnhatfr, abi_mean_fftr
 !use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,ider,idir,ipert,izero,my_natom,natom,nfft
 integer,intent(in)  :: usewvl
 integer,intent(in) :: nhatgrdim,nspden,ntypat
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: compch_fft
 type(pawang_type),intent(in) :: pawang

 integer,optional,intent(in) :: me_g0,mpi_comm_atom,mpi_comm_fft,mpi_comm_wvl,paral_kgb
 type(distribfft_type),optional,intent(in),target :: distribfft
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),qphon(3),rprimd(3,3),xred(3,natom)
 real(dp),intent(out) :: pawgrnhat(cplex*nfft,nspden,3*nhatgrdim)
 real(dp),intent(inout) :: pawnhat(cplex*nfft,nspden) !vz_i
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom),pawrhoij0(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

 integer,optional,target,intent(in) :: mpi_atmtab(:)

!Local variables ---------------------------------------
!scalars
 integer :: dplex,iatom,iatom_tot,ic,ierr,ii,ils,ilslm,irhoij,ispden,itypat
 integer :: jc,jrhoij,kc,klm,klmn,lmax,lmin,lm_size,mfgd,mm,mpi_comm_sphgrid
 integer :: my_comm_atom,my_comm_fft,nfgd,nfftot,option,optgr0,optgr1,optgr2,paral_kgb_fft
 logical :: compute_grad,compute_nhat,my_atmtab_allocated,need_frozen,paral_atom,qeq0
 logical :: compute_phonons,has_phase
 type(distribfft_type),pointer :: my_distribfft
 character(len = 500) :: msg
! type(mpi_type) :: mpi_enreg_fft
!arrays
 integer,pointer :: my_atmtab(:),my_atmtab_(:)
 real(dp) :: ro(cplex),ro_ql(cplex),tmp_compch_fft(nspden),tsec(2)
 real(dp),allocatable :: pawgrnhat_atm(:,:),pawnhat_atm(:),work(:,:)

! *************************************************************************

 ! DBG_ENTER("COLL")

 compute_nhat=(ider==0.or.ider==2)
 compute_grad=(ider==1.or.ider==2)
 compute_phonons=(ipert>0.and.ipert<=natom)

!Compatibility tests
 qeq0=(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15)
 if (present(mpi_comm_fft)) then
   if ((.not.present(paral_kgb)).or.(.not.present(me_g0))) then
      write(msg, "(A)") 'Need paral_kgb and me_g0 with mpi_comm_fft !'
     MSG_BUG(msg)
   end if
 end if
 if(ider>0.and.nhatgrdim==0) then
    write(msg, "(A)") ' Gradients of nhat required but not allocated !'
   MSG_BUG(msg)
 end if
 if (my_natom>0) then
   if(nspden>1.and.nspden/=pawrhoij(1)%nspden) then
      write(msg, "(A)") ' Wrong values for nspden and pawrhoij%nspden !'
     MSG_BUG(msg)
   end if
   if(nspden>1.and.nspden/=pawfgrtab(1)%nspden) then
      write(msg, "(A)") ' Wrong values for nspden and pawfgrtab%nspden !'
     MSG_BUG(msg)
   end if
   if(pawrhoij(1)%cplex<cplex) then
      write(msg, "(A)") '  Must have pawrhoij()%cplex >= cplex !'
     MSG_BUG(msg)
   end if
   if (compute_phonons.and.(.not.qeq0)) then
     if (pawfgrtab(1)%rfgd_allocated==0) then
        write(msg, "(A)") ' pawfgrtab()%rfgd array must be allocated  !'
       MSG_BUG(msg)
     end if
     if (compute_grad.and.(.not.compute_nhat)) then
        write(msg, "(A)") ' When q<>0, nhat gradients need nhat !'
       MSG_BUG(msg)
     end if
   end if
 end if

!nhat1 does not have to be computed for ddk or magn. field perturbation
 if (ipert==natom+1.or.ipert==natom+5) return

!Set up parallelism over atoms
 paral_atom=(present(mpi_comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xpaw_mpi_comm_self;if (present(mpi_comm_atom)) my_comm_atom=mpi_comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Initialisations
 if ((.not.compute_nhat).and.(.not.compute_grad)) return
 mfgd=zero;if (my_natom>0) mfgd=maxval(pawfgrtab(1:my_natom)%nfgd)
 dplex=cplex-1
 if (compute_nhat) then
   LIBPAW_ALLOCATE(pawnhat_atm,(cplex*mfgd))
   pawnhat=zero
 end if
 if (compute_grad) then
   LIBPAW_ALLOCATE(pawgrnhat_atm,(cplex*mfgd,3))
   pawgrnhat=zero
 end if

!------------------------------------------------------------------------
!----- Loop over atoms
!------------------------------------------------------------------------

 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

   itypat=pawrhoij(iatom)%itypat
   lm_size=pawfgrtab(iatom)%l_size**2
   need_frozen=((compute_nhat).and.(ipert==iatom_tot.or.ipert==natom+3.or.ipert==natom+4))
   nfgd=pawfgrtab(iatom)%nfgd
!  Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
   if (((compute_nhat).and.(pawfgrtab(iatom)%gylm_allocated==0)).or.&
&   ((compute_grad).and.(pawfgrtab(iatom)%gylmgr_allocated==0)).or.&
&   ((compute_grad.and.need_frozen).and.(pawfgrtab(iatom)%gylmgr2_allocated==0))) then
     optgr0=0;optgr1=0;optgr2=0
     if ((compute_nhat).and.(pawfgrtab(iatom)%gylm_allocated==0)) then
       if (allocated(pawfgrtab(iatom)%gylm))  then
         LIBPAW_DEALLOCATE(pawfgrtab(iatom)%gylm)
       end if
       LIBPAW_ALLOCATE(pawfgrtab(iatom)%gylm,(nfgd,pawfgrtab(iatom)%l_size**2))
       pawfgrtab(iatom)%gylm_allocated=2;optgr0=1
     end if
     if ((compute_grad).and.(pawfgrtab(iatom)%gylmgr_allocated==0)) then
       if (allocated(pawfgrtab(iatom)%gylmgr))  then
         LIBPAW_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
       end if
       LIBPAW_ALLOCATE(pawfgrtab(iatom)%gylmgr,(3,nfgd,pawfgrtab(iatom)%l_size**2))
       pawfgrtab(iatom)%gylmgr_allocated=2;optgr1=1
     end if
     if ((compute_grad.and.need_frozen).and.(pawfgrtab(iatom)%gylmgr2_allocated==0)) then
       if (allocated(pawfgrtab(iatom)%gylmgr2))  then
         LIBPAW_DEALLOCATE(pawfgrtab(iatom)%gylmgr2)
       end if
       LIBPAW_ALLOCATE(pawfgrtab(iatom)%gylmgr2,(6,nfgd,pawfgrtab(iatom)%l_size**2))
       pawfgrtab(iatom)%gylmgr2_allocated=2;optgr2=1
     end if
     if (optgr0+optgr1+optgr2>0) then
       call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,pawfgrtab(iatom)%gylmgr2,&
&       lm_size,nfgd,optgr0,optgr1,optgr2,pawtab(itypat),pawfgrtab(iatom)%rfgd)
     end if
   end if


!  Eventually compute exp(-i.q.r) factors for the current atom (if not already done)
   if (compute_phonons.and.(.not.qeq0).and.pawfgrtab(iatom)%expiqr_allocated==0) then
     if (allocated(pawfgrtab(iatom)%expiqr))  then
       LIBPAW_DEALLOCATE(pawfgrtab(iatom)%expiqr)
     end if
     LIBPAW_ALLOCATE(pawfgrtab(iatom)%expiqr,(2,nfgd))
     call pawexpiqr(pawfgrtab(iatom)%expiqr,gprimd,nfgd,qphon,&
&     pawfgrtab(iatom)%rfgd,xred(:,iatom_tot))
     pawfgrtab(iatom)%expiqr_allocated=2
   end if
   has_phase=(compute_phonons.and.pawfgrtab(iatom)%expiqr_allocated/=0)

!  Eventually compute frozen part of nhat for the current atom (if not already done)
   if ((need_frozen).and.((pawfgrtab(iatom)%nhatfr_allocated==0).or.&
&   (compute_grad.and.pawfgrtab(iatom)%nhatfrgr_allocated==0))) then
     if (allocated(pawfgrtab(iatom)%nhatfr))  then
       LIBPAW_DEALLOCATE(pawfgrtab(iatom)%nhatfr)
     end if
     LIBPAW_ALLOCATE(pawfgrtab(iatom)%nhatfr,(nfgd,pawfgrtab(iatom)%nspden))
     option=0;pawfgrtab(iatom)%nhatfr_allocated=2
     if (compute_grad) then
       option=1
       if (allocated(pawfgrtab(iatom)%nhatfrgr))  then
         LIBPAW_DEALLOCATE(pawfgrtab(iatom)%nhatfrgr)
       end if
       LIBPAW_ALLOCATE(pawfgrtab(iatom)%nhatfrgr,(3,nfgd,pawfgrtab(iatom)%nspden))
       pawfgrtab(iatom)%nhatfrgr_allocated=2
     end if
     nullify(my_atmtab_);if (associated(my_atmtab)) my_atmtab_ => my_atmtab
     call abi_pawnhatfr(option,idir,ipert,my_natom,natom,nspden,ntypat,pawang,pawfgrtab,&
&     pawrhoij0,pawtab,rprimd,mpi_comm_atom=my_comm_atom,mpi_atmtab=my_atmtab_)
   end if

!  ------------------------------------------------------------------------
!  ----- Loop over density components
!  ------------------------------------------------------------------------

   do ispden=1,nspden

     if (compute_nhat) pawnhat_atm(1:cplex*nfgd)=zero
     if (compute_grad) pawgrnhat_atm(1:cplex*nfgd,1:3)=zero

!    ------------------------------------------------------------------------
!    ----- Loop over ij channels (basis components)
!    ------------------------------------------------------------------------
     jrhoij=1
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
       klmn=pawrhoij(iatom)%rhoijselect(irhoij)
       klm =pawtab(itypat)%indklmn(1,klmn)
       lmin=pawtab(itypat)%indklmn(3,klmn)
       lmax=pawtab(itypat)%indklmn(4,klmn)

!      Retrieve rhoij
       if (nspden/=2) then
         ro(1:cplex)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
       else
         if (ispden==1) then
           ro(1:cplex)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,1)&
&           +pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,2)
         else if (ispden==2) then
           ro(1:cplex)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,1)
         end if
       end if
       ro(1:cplex)=pawtab(itypat)%dltij(klmn)*ro(1:cplex)
       if (compute_nhat) then
         if (cplex==1) then
           do ils=lmin,lmax,2
             do mm=-ils,ils
               ilslm=ils*ils+ils+mm+1
               if (pawang%gntselect(ilslm,klm)>0) then
                 ro_ql(1)=ro(1)*pawtab(itypat)%qijl(ilslm,klmn)
                 do ic=1,nfgd
                   pawnhat_atm(ic)=pawnhat_atm(ic)+ro_ql(1)*pawfgrtab(iatom)%gylm(ic,ilslm)
                 end do
               end if
             end do
           end do
         else
           do ils=lmin,lmax,2
             do mm=-ils,ils
               ilslm=ils*ils+ils+mm+1
               if (pawang%gntselect(ilslm,klm)>0) then
                 ro_ql(1:2)=ro(1:2)*pawtab(itypat)%qijl(ilslm,klmn)
                 do ic=1,nfgd
                   jc=2*ic-1
                   pawnhat_atm(jc:jc+1)=pawnhat_atm(jc:jc+1)+ro_ql(1:2)*pawfgrtab(iatom)%gylm(ic,ilslm)
                 end do
               end if
             end do
           end do
         end if
       end if

       if (compute_grad) then
         if (cplex==1) then
           do ils=lmin,lmax,2
             do mm=-ils,ils
               ilslm=ils*ils+ils+mm+1
               if (pawang%gntselect(ilslm,klm)>0) then
                 ro_ql(1)=ro(1)*pawtab(itypat)%qijl(ilslm,klmn)
                 do ic=1,nfgd
                   pawgrnhat_atm(ic,1)=pawgrnhat_atm(ic,1)+ro_ql(1)*pawfgrtab(iatom)%gylmgr(1,ic,ilslm)
                   pawgrnhat_atm(ic,2)=pawgrnhat_atm(ic,2)+ro_ql(1)*pawfgrtab(iatom)%gylmgr(2,ic,ilslm)
                   pawgrnhat_atm(ic,3)=pawgrnhat_atm(ic,3)+ro_ql(1)*pawfgrtab(iatom)%gylmgr(3,ic,ilslm)
                 end do
               end if
             end do
           end do
         else
           do ils=lmin,lmax,2
             do mm=-ils,ils
               ilslm=ils*ils+ils+mm+1
               if (pawang%gntselect(ilslm,klm)>0) then
                 ro_ql(1:2)=ro(1:2)*pawtab(itypat)%qijl(ilslm,klmn)
                 do ic=1,nfgd
                   jc=2*ic-1
                   pawgrnhat_atm(jc:jc+1,1)=pawgrnhat_atm(jc:jc+1,1) &
&                   +ro_ql(1:2)*pawfgrtab(iatom)%gylmgr(1,ic,ilslm)
                   pawgrnhat_atm(jc:jc+1,2)=pawgrnhat_atm(jc:jc+1,2) &
&                   +ro_ql(1:2)*pawfgrtab(iatom)%gylmgr(2,ic,ilslm)
                   pawgrnhat_atm(jc:jc+1,3)=pawgrnhat_atm(jc:jc+1,3) &
&                   +ro_ql(1:2)*pawfgrtab(iatom)%gylmgr(3,ic,ilslm)
                 end do
               end if
             end do
           end do
         end if
       end if

!      ------------------------------------------------------------------------
!      ----- End loop over ij channels
!      ------------------------------------------------------------------------
       jrhoij=jrhoij+pawrhoij(iatom)%cplex
     end do

!    If RF calculation, add frozen part of 1st-order compensation density
     if (need_frozen) then
       if (cplex==1) then
         do ic=1,nfgd
           pawnhat_atm(ic)=pawnhat_atm(ic)+pawfgrtab(iatom)%nhatfr(ic,ispden)
         end do
       else
         do ic=1,nfgd
           jc=2*ic-1
           pawnhat_atm(jc)=pawnhat_atm(jc)+pawfgrtab(iatom)%nhatfr(ic,ispden)
         end do
       end if
       if (compute_grad) then
         if (cplex==1) then
           do ic=1,nfgd
             pawgrnhat_atm(ic,1)=pawgrnhat_atm(ic,1)+pawfgrtab(iatom)%nhatfrgr(1,ic,ispden)
             pawgrnhat_atm(ic,2)=pawgrnhat_atm(ic,2)+pawfgrtab(iatom)%nhatfrgr(2,ic,ispden)
             pawgrnhat_atm(ic,3)=pawgrnhat_atm(ic,3)+pawfgrtab(iatom)%nhatfrgr(3,ic,ispden)
           end do
         else
           do ic=1,nfgd
             jc=2*ic-1
             pawgrnhat_atm(jc,1)=pawgrnhat_atm(jc,1)+pawfgrtab(iatom)%nhatfrgr(1,ic,ispden)
             pawgrnhat_atm(jc,2)=pawgrnhat_atm(jc,2)+pawfgrtab(iatom)%nhatfrgr(2,ic,ispden)
             pawgrnhat_atm(jc,3)=pawgrnhat_atm(jc,3)+pawfgrtab(iatom)%nhatfrgr(3,ic,ispden)
           end do
         end if
       end if
     end if

!    If needed, multiply eventually by exp(-i.q.r) phase
     if (has_phase) then
       if (cplex==1) then
         do ic=1,nfgd
           pawnhat_atm(ic)=pawnhat_atm(ic)*pawfgrtab(iatom)%expiqr(1,ic)
         end do
       else
         do ic=1,nfgd
           jc=2*ic-1
           ro_ql(1)= pawfgrtab(iatom)%expiqr(1,ic)
           ro_ql(2)=-pawfgrtab(iatom)%expiqr(2,ic)
           ro(1:2)=pawnhat_atm(jc:jc+1)
           pawnhat_atm(jc  )=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
           pawnhat_atm(jc+1)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
         end do
       end if
       if (compute_grad) then
         if (cplex==1) then
           do ic=1,nfgd
             pawgrnhat_atm(ic,1:3)=pawgrnhat_atm(ic,1:3)*pawfgrtab(iatom)%expiqr(1,ic)
           end do
         else
           do ic=1,nfgd
             jc=2*ic-1
!            dn^hat(r)/dr_i * exp(-i.q.r)
             ro_ql(1)= pawfgrtab(iatom)%expiqr(1,ic)
             ro_ql(2)=-pawfgrtab(iatom)%expiqr(2,ic)
             do ii=1,3
               ro(1:2)=pawgrnhat_atm(jc:jc+1,ii)
               pawgrnhat_atm(jc  ,ii)=ro(1)*ro_ql(1)-ro(2)*ro_ql(2)
               pawgrnhat_atm(jc+1,ii)=ro(2)*ro_ql(1)+ro(1)*ro_ql(2)
             end do
!            -i.q_i * [n^hat(r).exp(-i.q.r)]
             ro(1:2)=pawnhat_atm(jc:jc+1)
             do ii=1,3
               pawgrnhat_atm(jc  ,ii)=pawgrnhat_atm(kc  ,ii)+qphon(ii)*ro(2)
               pawgrnhat_atm(jc+1,ii)=pawgrnhat_atm(kc+1,ii)-qphon(ii)*ro(1)
             end do
           end do
         end if
       end if
     end if

!    Add the contribution of the atom to the compensation charge
     if (compute_nhat) then
       if (cplex==1) then
         do ic=1,nfgd
           kc=pawfgrtab(iatom)%ifftsph(ic)
           pawnhat(kc,ispden)=pawnhat(kc,ispden)+pawnhat_atm(ic)
         end do
       else
         do ic=1,nfgd
           jc=2*ic-1;kc=2*pawfgrtab(iatom)%ifftsph(ic)-1
           pawnhat(kc:kc+1,ispden)=pawnhat(kc:kc+1,ispden)+pawnhat_atm(jc:jc+1)
         end do
       end if
     end if
     if (compute_grad) then
       if (cplex==1) then
         do ic=1,nfgd
           kc=pawfgrtab(iatom)%ifftsph(ic)
           pawgrnhat(kc,ispden,1:3)=pawgrnhat(kc,ispden,1:3)+pawgrnhat_atm(ic,1:3)
         end do
       else
         do ic=1,nfgd
           jc=2*ic-1;kc=2*pawfgrtab(iatom)%ifftsph(ic)-1
           do ii=1,3
             pawgrnhat(kc:kc+1,ispden,ii)=pawgrnhat(kc:kc+1,ispden,ii)+pawgrnhat_atm(jc:jc+1,ii)
           end do
         end do
       end if
     end if
!    ------------------------------------------------------------------------
!    ----- End loop over density components
!    ------------------------------------------------------------------------
   end do

   if (pawfgrtab(iatom)%gylm_allocated==2) then
     LIBPAW_DEALLOCATE(pawfgrtab(iatom)%gylm)
     LIBPAW_ALLOCATE(pawfgrtab(iatom)%gylm,(0,0))
     pawfgrtab(iatom)%gylm_allocated=0
   end if
   if (pawfgrtab(iatom)%gylmgr_allocated==2) then
     LIBPAW_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
     LIBPAW_ALLOCATE(pawfgrtab(iatom)%gylmgr,(0,0,0))
     pawfgrtab(iatom)%gylmgr_allocated=0
   end if
   if (pawfgrtab(iatom)%gylmgr2_allocated==2) then
     LIBPAW_DEALLOCATE(pawfgrtab(iatom)%gylmgr2)
     LIBPAW_ALLOCATE(pawfgrtab(iatom)%gylmgr2,(0,0,0))
     pawfgrtab(iatom)%gylmgr2_allocated=0
   end if
   if (pawfgrtab(iatom)%nhatfr_allocated==2) then
     LIBPAW_DEALLOCATE(pawfgrtab(iatom)%nhatfr)
     LIBPAW_ALLOCATE(pawfgrtab(iatom)%nhatfr,(0,0))
     pawfgrtab(iatom)%nhatfr_allocated=0
   end if
   if (pawfgrtab(iatom)%nhatfrgr_allocated==2) then
     LIBPAW_DEALLOCATE(pawfgrtab(iatom)%nhatfrgr)
     LIBPAW_ALLOCATE(pawfgrtab(iatom)%nhatfrgr,(0,0,0))
     pawfgrtab(iatom)%nhatfrgr_allocated=0
   end if
   if (pawfgrtab(iatom)%expiqr_allocated==2) then
     LIBPAW_DEALLOCATE(pawfgrtab(iatom)%expiqr)
     LIBPAW_ALLOCATE(pawfgrtab(iatom)%expiqr,(0,0))
     pawfgrtab(iatom)%expiqr_allocated=0
   end if

!  ------------------------------------------------------------------------
!  ----- End loop over atoms
!  ------------------------------------------------------------------------
 end do

!----- Free some memory
 if (compute_nhat) then
   LIBPAW_DEALLOCATE(pawnhat_atm)
 end if
 if (compute_grad) then
   LIBPAW_DEALLOCATE(pawgrnhat_atm)
 end if

!----- Reduction in case of parallelism
 if (paral_atom) then
   !call timab(48,1,tsec)
   if (compute_nhat) then
     call xpaw_mpi_sum(pawnhat,my_comm_atom,ierr)
   end if
   if (compute_grad) then
     call xpaw_mpi_sum(pawgrnhat,my_comm_atom,ierr)
   end if
   !call timab(48,2,tsec)
 end if

!----- Avoid unbalanced g-components numerical errors
!!$ if (izero==1.and.compute_nhat.and.usewvl==0) then
!!$!  Create fake mpi_enreg to wrap fourdp
!!$   if (present(distribfft)) then
!!$     my_distribfft => distribfft
!!$   else
!!$     LIBPAW_DATATYPE_ALLOCATE(my_distribfft,)
!!$     call init_distribfft_seq(my_distribfft,'f',ngfft(2),ngfft(3),'fourdp')
!!$   end if
!!$   call initmpi_seq(mpi_enreg_fft)
!!$   LIBPAW_DATATYPE_DEALLOCATE(mpi_enreg_fft%distribfft)
!!$   if (present(mpi_comm_fft)) then
!!$     call set_mpi_enreg_fft(mpi_enreg_fft,mpi_comm_fft,my_distribfft,me_g0,paral_kgb)
!!$     my_comm_fft=mpi_comm_fft;paral_kgb_fft=paral_kgb
!!$   else
!!$     my_comm_fft=xpaw_mpi_comm_self;paral_kgb_fft=0;
!!$     mpi_enreg_fft%distribfft => my_distribfft
!!$   end if
!!$!  do FFT
!!$   LIBPAW_ALLOCATE(work,(2,nfft))
!!$   do ispden=1,min(2,nspden)
!!$     call fourdp(cplex,work,pawnhat(:,ispden),-1,mpi_enreg_fft,nfft,ngfft,paral_kgb_fft,0)
!!$     call zerosym(work,2,ngfft(1),ngfft(2),ngfft(3),mpi_comm_fft=my_comm_fft,distribfft=my_distribfft)
!!$     call fourdp(cplex,work,pawnhat(:,ispden),+1,mpi_enreg_fft,nfft,ngfft,paral_kgb_fft,0)
!!$   end do
!!$   LIBPAW_DEALLOCATE(work)
!!$!  Destroy fake mpi_enreg
!!$   call unset_mpi_enreg_fft(mpi_enreg_fft)
!!$   if (.not.present(distribfft)) then
!!$     call destroy_distribfft(my_distribfft)
!!$     LIBPAW_DATATYPE_DEALLOCATE(my_distribfft)
!!$   end if
!!$ end if

!----- Computation of compensation charge over real space grid
 if (compute_nhat.and.ipert==0) then
    !mpi communicators for spherical grid:
    mpi_comm_sphgrid=xpaw_mpi_comm_self !no communicators passed
    if(present(mpi_comm_fft) .and. usewvl==0) mpi_comm_sphgrid=mpi_comm_fft
    if(present(mpi_comm_wvl) .and. usewvl==1) mpi_comm_sphgrid=mpi_comm_wvl
    
    nfftot=PRODUCT(ngfft(1:3))
    call abi_mean_fftr(pawnhat,tmp_compch_fft,nfft,nfftot,1,mpi_comm_sphgrid)
    compch_fft = tmp_compch_fft(1)
    compch_fft=compch_fft*ucvol
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 ! DBG_EXIT("COLL")

end subroutine abi_pawmknhat
!!***

!!****f* m_cgtools/abi_mean_fftr
!! NAME
!! abi_mean_fftr
!!
!! FUNCTION
!!  Compute the mean of an arraysp(nfft,nspden), over the FFT grid, for each component nspden, 
!!  and return it in meansp(nspden). 
!!  Take into account the spread of the array due to parallelism: the actual number of fft
!!  points is nfftot, but the number of points on this proc is nfft only.
!!  So : for ispden from 1 to nspden
!!       meansp(ispden) = sum(ifft=1,nfftot) arraysp(ifft,ispden) / nfftot
!!
!! INPUTS
!!  arraysp(nfft,nspden)=the array whose average has to be computed
!!  nfft=number of FFT points stored by one proc
!!  nfftot=total number of FFT points
!!  nspden=number of spin-density components
!!
!! OUTPUT
!!  meansp(nspden)=mean value for each nspden component
!!
!! PARENTS
!!      fresid,multipoles_fftr,newvtr,abi_pawmknhat,prcref,prcref_PMA
!!      psolver_rhohxc,rhohxc,rhohxcpositron,rhotov
!!
!! CHILDREN
!!      timab,xmpi_sum
!!
!! SOURCE

subroutine abi_mean_fftr(arraysp,meansp,nfft,nfftot,nspden,mpi_comm_sphgrid)

 USE_DEFS
 USE_MPI_WRAPPERS

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_mean_fftr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nfftot,nspden
 integer,intent(in),optional:: mpi_comm_sphgrid
!arrays
 real(dp),intent(in) :: arraysp(nfft,nspden)
 real(dp),intent(out) :: meansp(nspden)

!Local variables-------------------------------
!scalars
 integer :: ierr,ifft,ispden,nproc_sphgrid
 real(dp) :: invnfftot,tmean

! *************************************************************************

 invnfftot=one/(dble(nfftot))

 do ispden=1,nspden
   tmean=zero
!$OMP PARALLEL DO REDUCTION(+:tmean)
   do ifft=1,nfft
     tmean=tmean+arraysp(ifft,ispden)
   end do
   meansp(ispden)=tmean*invnfftot
 end do

!XG030514 : MPIWF The values of meansp(ispden) should
!now be summed accross processors in the same WF group, and spread on all procs.
!if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_fft; reduce=.true.
 if(present(mpi_comm_sphgrid)) then
   nproc_sphgrid=xpaw_mpi_comm_size(mpi_comm_sphgrid)
   if(nproc_sphgrid>1) then
     call xpaw_mpi_sum(meansp,mpi_comm_sphgrid,ierr)
   end if
 end if

end subroutine abi_mean_fftr
!!***
