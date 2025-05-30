!> @file
!!  Time-Dependent DFT ai la Casida
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Calculate the coupling matrix for the TD-DFT a la Casida
subroutine tddft_casida(iproc,nproc,atoms,rxyz,hxh,hyh,hzh,n3p,n3parr,Glr,tddft_approach,orbs,&
     orbsv,i3s,fxc,pkernelseq,psi,psiv,exc_fac,bitp)
  use module_base
  use module_types
  use locregs
  use box
  implicit none
  integer, intent(in) :: iproc,nproc,n3p,i3s
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: atoms
  type(orbitals_data), intent(in) :: orbs,orbsv
  type(locreg_descriptors), intent(in) :: Glr
  character(len=4), intent(in) :: tddft_approach
  integer, dimension(0:nproc-1), intent(in) :: n3parr
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  real(wp), dimension(Glr%d%n1i,Glr%d%n2i,n3p), intent(in) :: fxc
  type(coulomb_operator), intent(inout) :: pkernelseq
  real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
  real(wp), dimension(orbsv%npsidim_orbs), intent(in) :: psiv
  real(gp), intent(in) :: exc_fac
  type(box_iterator) :: bitp
  !local variables
  real(gp), dimension(3) :: chargec
  real(wp), dimension(:), allocatable :: psirocc,psirvirt

  !temporary call to the coupling matrix calculation
  psirocc = f_malloc(max(max(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*orbs%norbp, n3parr(0)*orbs%norb), 1),id='psirocc')
  psirvirt = f_malloc(max(max(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*orbsv%norbp, n3parr(0)*orbsv%norb), 1),id='psirvirt')

  call prepare_psirocc(iproc,nproc,Glr,orbs,n3p,n3parr,psi,psirocc)

  call prepare_psirocc(iproc,nproc,Glr,orbsv,n3p,n3parr,psiv,psirvirt)

  call center_of_charge(atoms,rxyz,chargec)

  call calculate_coupling_matrix(iproc,nproc,bitp,tddft_approach,orbs%nspin,Glr%d%n1i*Glr%d%n2i*n3p,orbs,orbsv,&
       chargec,pkernelseq,fxc,psirocc,psirvirt)

!!$  call coupling_matrix_prelim(iproc,nproc,atoms%astruct%geocode,tddft_approach,orbs%nspin,Glr,orbs,orbsv,&
!!$       i3s,n3p,hxh,hyh,hzh,chargec,pkernelseq,fxc,psirocc,psirvirt,exc_fac)

  call f_free(psirocc)

  call f_free(psirvirt)

end subroutine tddft_casida
