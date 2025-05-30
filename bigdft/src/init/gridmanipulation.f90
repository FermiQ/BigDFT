!> @file
!!  Routines to manipulate the grid
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calculates the overall size of the simulation cell 
!! and shifts the atoms such that their position is the most symmetric possible.
!! Assign these values to the global localisation region descriptor.
subroutine system_size(atoms,rxyz,crmult,frmult,hx,hy,hz,OCLconv,Glr)
   use module_base
   use module_types
   use yaml_strings, only: yaml_toa
   use locregs
   implicit none
   type(atoms_data), intent(inout) :: atoms
   real(gp), intent(in) :: crmult,frmult
   real(gp), dimension(3,atoms%astruct%nat), intent(inout) :: rxyz
   !real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
   real(gp), intent(inout) :: hx,hy,hz
   logical, intent(in) :: OCLconv
   type(locreg_descriptors), intent(out) :: Glr
   !Local variables
   !character(len=*), parameter :: subname='system_size'
   integer, parameter :: lupfil=14
   real(gp), parameter :: eps_mach=1.e-12_gp
   integer :: iat,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3!,n1i,n2i,n3i
   real(gp) :: ri,rad,cxmin,cxmax,cymin,cymax,czmin,czmax,alatrue1,alatrue2,alatrue3
   real(gp), dimension(3) :: hgridsh

   !check the geometry code with the grid spacings
   if (atoms%astruct%geocode == 'F' .and. (hx/=hy .or. hx/=hz .or. hy/=hz)) then
      call f_err_throw('Grid spacings must be equal' // &
           ' in the Free BC case, while hgrids = '//&
           trim(yaml_toa((/ hx, hy, hz/),fmt='(f7.4)')),&
           err_name='BIGDFT_INPUT_VARIABLES_ERROR')
      !write(*,'(1x,a,3(1x,F7.4))') 'ERROR: The values of the grid spacings must be equal' // &
      !     & ' in the Free BC case, while hgrids = ', hx, hy, hz
      return
   end if

   !Special case if no atoms (no posinp by error or electron gas)
   if (atoms%astruct%nat == 0) then
      ri = 0.0_gp
   else
      ri = 1.e10_gp
   end if

   !calculate the extremes of the boxes taking into account the spheres around the atoms
   cxmax = -ri 
   cxmin =  ri

   cymax = -ri 
   cymin =  ri

   czmax = -ri 
   czmin =  ri

   do iat=1,atoms%astruct%nat

      rad=atoms%radii_cf(atoms%astruct%iatype(iat),1)*crmult

      cxmax=max(cxmax,rxyz(1,iat)+rad) 
      cxmin=min(cxmin,rxyz(1,iat)-rad)

      cymax=max(cymax,rxyz(2,iat)+rad) 
      cymin=min(cymin,rxyz(2,iat)-rad)

      czmax=max(czmax,rxyz(3,iat)+rad) 
      czmin=min(czmin,rxyz(3,iat)-rad)
   enddo

   !eliminate epsilon form the grid size calculation
   !!  cxmax=cxmax+eps_mach 
   !!  cymax=cymax+eps_mach  
   !!  czmax=czmax+eps_mach  
   !!
   !!  cxmin=cxmin-eps_mach
   !!  cymin=cymin-eps_mach
   !!  czmin=czmin-eps_mach

   !define the box sizes for free BC, and calculate dimensions for the fine grid with ISF
   select case (atoms%astruct%geocode)
   
   case('F')
      atoms%astruct%cell_dim(1)=(cxmax-cxmin)
      atoms%astruct%cell_dim(2)=(cymax-cymin)
      atoms%astruct%cell_dim(3)=(czmax-czmin)

      ! grid sizes n1,n2,n3
      n1=int(atoms%astruct%cell_dim(1)/hx)
      !if (mod(n1,2)==1) n1=n1+1
      n2=int(atoms%astruct%cell_dim(2)/hy)
      !if (mod(n2,2)==1) n2=n2+1
      n3=int(atoms%astruct%cell_dim(3)/hz)
      !if (mod(n3,2)==1) n3=n3+1
      alatrue1=real(n1,gp)*hx
      alatrue2=real(n2,gp)*hy
      alatrue3=real(n3,gp)*hz

!!$      n1i=2*n1+31
!!$      n2i=2*n2+31
!!$      n3i=2*n3+31

   case('P')
      !define the grid spacings, controlling the FFT compatibility
      call correct_grid(atoms%astruct%cell_dim(1),hx,n1)
      call correct_grid(atoms%astruct%cell_dim(2),hy,n2)
      call correct_grid(atoms%astruct%cell_dim(3),hz,n3)
      alatrue1=(cxmax-cxmin)
      alatrue2=(cymax-cymin)
      alatrue3=(czmax-czmin)

!!$      n1i=2*n1+2
!!$      n2i=2*n2+2
!!$      n3i=2*n3+2

   case('S')
      call correct_grid(atoms%astruct%cell_dim(1),hx,n1)
      atoms%astruct%cell_dim(2)=(cymax-cymin)
      call correct_grid(atoms%astruct%cell_dim(3),hz,n3)

      alatrue1=(cxmax-cxmin)
      n2=int(atoms%astruct%cell_dim(2)/hy)
      alatrue2=real(n2,gp)*hy
      alatrue3=(czmax-czmin)

!!$      n1i=2*n1+2
!!$      n2i=2*n2+31
!!$      n3i=2*n3+2

   case default
      call f_err_throw('Illegal geocode in system_size',err_id=BIGDFT_INPUT_VARIABLES_ERROR)

   end select

   !balanced shift taking into account the missing space
   cxmin=cxmin+0.5_gp*(atoms%astruct%cell_dim(1)-alatrue1)
   cymin=cymin+0.5_gp*(atoms%astruct%cell_dim(2)-alatrue2)
   czmin=czmin+0.5_gp*(atoms%astruct%cell_dim(3)-alatrue3)

   !correct the box sizes for the isolated case
   select case(atoms%astruct%geocode)
   case('F')
      atoms%astruct%cell_dim(1)=alatrue1
      atoms%astruct%cell_dim(2)=alatrue2
      atoms%astruct%cell_dim(3)=alatrue3
   case('S')
      cxmin=0.0_gp
      atoms%astruct%cell_dim(2)=alatrue2
      czmin=0.0_gp
   case('P')
      !for the moment we do not put the shift, at the end it will be tested
      !here we should put the center of mass
      cxmin=0.0_gp
      cymin=0.0_gp
      czmin=0.0_gp
   end select

   !assign the shift to the atomic positions
   atoms%astruct%shift(1)=cxmin
   atoms%astruct%shift(2)=cymin
   atoms%astruct%shift(3)=czmin

   !here we can put a modulo operation for periodic directions
   do iat=1,atoms%astruct%nat
      rxyz(1,iat)=rxyz(1,iat)-atoms%astruct%shift(1)
      rxyz(2,iat)=rxyz(2,iat)-atoms%astruct%shift(2)
      rxyz(3,iat)=rxyz(3,iat)-atoms%astruct%shift(3)
   enddo

   ! fine grid size (needed for creation of input wavefunction, preconditioning)
   if (atoms%astruct%nat == 0) then
      !For homogeneous gaz, we fill the box with the fine grid
      nfl1=0 
      nfl2=0 
      nfl3=0

      nfu1=n1
      nfu2=n2
      nfu3=n3
   else
      !we start with nfl max to find th emin and nfu min to find the max
      nfl1=n1 
      nfl2=n2 
      nfl3=n3

      nfu1=0 
      nfu2=0 
      nfu3=0
   end if

   do iat=1,atoms%astruct%nat
      rad=atoms%radii_cf(atoms%astruct%iatype(iat),2)*frmult
      if (rad > 0.0_gp) then
         nfl1=min(nfl1,ceiling((rxyz(1,iat)-rad)/hx - eps_mach))
         nfu1=max(nfu1,floor((rxyz(1,iat)+rad)/hx + eps_mach))

         nfl2=min(nfl2,ceiling((rxyz(2,iat)-rad)/hy - eps_mach))
         nfu2=max(nfu2,floor((rxyz(2,iat)+rad)/hy + eps_mach))

         nfl3=min(nfl3,ceiling((rxyz(3,iat)-rad)/hz - eps_mach)) 
         nfu3=max(nfu3,floor((rxyz(3,iat)+rad)/hz + eps_mach))
      end if
   enddo

   !correct the values of the delimiter if they go outside the box
   if (nfl1 < 0 .or. nfu1 > n1) then
      nfl1=0
      nfu1=n1
   end if
   if (nfl2 < 0 .or. nfu2 > n2) then
      nfl2=0
      nfu2=n2
   end if
   if (nfl3 < 0 .or. nfu3 > n3) then
      nfl3=0
      nfu3=n3
   end if

   !correct the values of the delimiter if there are no wavelets
   if (nfl1 == n1 .and. nfu1 == 0) then
      nfl1=n1/2
      nfu1=n1/2
   end if
   if (nfl2 == n2 .and. nfu2 == 0) then
      nfl2=n2/2
      nfu2=n2/2
   end if
   if (nfl3 == n3 .and. nfu3 == 0) then
      nfl3=n3/2
      nfu3=n3/2
   end if

   hgridsh(1)=0.5_gp*hx
   hgridsh(2)=0.5_gp*hy
   hgridsh(3)=0.5_gp*hz

   !assign the values
   call init_lr(Glr,atoms%astruct%geocode,hgridsh,n1,n2,n3,&
        nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,&
        hybrid_flag=.not. OCLconv)

END SUBROUTINE system_size


!> Here the dimensions should be corrected in order to 
!! allow the fft for the preconditioner and for Poisson Solver
subroutine correct_grid(a,h,n)
   use module_base
   use Poisson_Solver, except_dp => dp, except_gp => gp
   implicit none
   real(gp), intent(in) :: a
   integer, intent(inout) :: n
   real(gp), intent(inout) :: h
   !local variables
   integer :: m,m2,nt

   n=ceiling(a/h)-1
   nt=n+1
   do
      !correct the direct dimension
      call fourier_dim(nt,m)

      !control if the double of this dimension is compatible with the FFT
      call fourier_dim(2*m,m2)
      !if this check is passed both the preconditioner and the PSolver works
      if (m2==2*m .and. mod(m,2) ==0) exit !only even dimensions are considered so far

      nt=m+1
   end do
   n=m-1

   !!!  !here the dimensions should be corrected in order to 
   !!!  !allow the fft for the preconditioner
   !!!  m=2*n+2
   !!!  do 
   !!!     call fourier_dim(m,m)
   !!!     if ((m/2)*2==m) then
   !!!        n=(m-2)/2
   !!!        exit
   !!!     else
   !!!        m=m+1
   !!!     end if
   !!!  end do

   h=a/real(n+1,gp)

END SUBROUTINE correct_grid


!> Calculates the length of the keys describing a wavefunction data structure
subroutine num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
   use dynamic_memory
   implicit none
   integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3
   logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid 
   integer, intent(out) :: mseg,mvctr
   !local variables
   logical :: plogrid
   integer :: i1,i2,i3,nsrt,nend,nsrti,nendi,mvctri

   call f_routine(id='num_segkeys')

   mvctr=0
   nsrt=0
   nend=0
   !$omp parallel default(private) shared(nl3,nu3,nl2,nu2,nl1,nu1,logrid,mvctr,nsrt,nend)
   mvctri=0
   nsrti=0
   nendi=0
   !$omp do  
   do i3=nl3,nu3 
      do i2=nl2,nu2
         plogrid=.false.
         do i1=nl1,nu1
            if (logrid(i1,i2,i3)) then
               mvctri=mvctri+1
               if (.not. plogrid) then
                  nsrti=nsrti+1
               endif
            else
               if (plogrid) then
                  nendi=nendi+1
               endif
            endif
            plogrid=logrid(i1,i2,i3)
         enddo
         if (plogrid) then
            nendi=nendi+1
         endif
      enddo
   enddo
   !$omp enddo
   !$omp critical
   mvctr=mvctr+mvctri
   nsrt=nsrt+nsrti
   nend=nend+nendi
   !$omp end critical
   !$omp end parallel
   if (nend /= nsrt) then 
      write(*,*)' ERROR: nend <> nsrt',nend,nsrt
      stop 
   endif
   mseg=nend

   call f_release_routine()

END SUBROUTINE num_segkeys


!> Calculates the keys describing a wavefunction data structure
subroutine segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,keyg,keyv)
   use dynamic_memory
   implicit none
   integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,mseg
   logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid  
   integer, dimension(mseg), intent(out) :: keyv
   integer, dimension(2,mseg), intent(out) :: keyg
   !local variables
   logical :: plogrid
   integer :: mvctr,nsrt,nend,i1,i2,i3,ngridp,np,n1p1

   call f_routine(id='segkeys')

   mvctr=0
   nsrt=0
   nend=0
   n1p1=n1+1
   np=n1p1*(n2+1)
   do i3=nl3,nu3 
      do i2=nl2,nu2
         plogrid=.false.
         do i1=nl1,nu1
            ngridp=i3*np + i2*n1p1 + i1+1
            if (logrid(i1,i2,i3)) then
               mvctr=mvctr+1
               if (.not. plogrid) then
                  nsrt=nsrt+1
                  keyg(1,nsrt)=ngridp
                  keyv(nsrt)=mvctr
               endif
            else
               if (plogrid) then
                  nend=nend+1
                  keyg(2,nend)=ngridp-1
               endif
            endif
            plogrid=logrid(i1,i2,i3)
         enddo
         if (plogrid) then
            nend=nend+1
            keyg(2,nend)=ngridp
         endif
      enddo
   enddo
   if (nend /= nsrt) then 
      write(*,*) 'nend , nsrt',nend,nsrt
      stop 'nend <> nsrt'
   endif
   !mseg=nend

   call f_release_routine()

END SUBROUTINE segkeys


subroutine export_grids(fname, atoms, rxyz, hx, hy, hz, n1, n2, n3, logrid_c, logrid_f)
  use module_defs, only: gp
  use module_types
  implicit none
  character(len = *), intent(in) :: fname
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
  real(gp), intent(in) :: hx, hy, hz
  integer, intent(in) :: n1, n2, n3
  logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid_c
  logical, dimension(0:n1,0:n2,0:n3), intent(in), optional :: logrid_f

  integer :: nvctr, iat, i3, i2, i1

  nvctr = 0
  do i3=0,n3  
     do i2=0,n2  
        do i1=0,n1
           if (logrid_c(i1,i2,i3)) nvctr = nvctr + 1
        enddo
     enddo
  end do
  if (present(logrid_f)) then
     do i3=0,n3  
        do i2=0,n2  
           do i1=0,n1
              if (logrid_f(i1,i2,i3)) nvctr = nvctr + 1
           enddo
        enddo
     end do
  end if

  ! Create the file grid.xyz to visualize the grid of functions
  open(unit=22,file=fname,status='unknown')
  write(22,*) nvctr+atoms%astruct%nat,' atomic'
  if (atoms%astruct%geocode=='F') then
     write(22,*)'complete simulation grid with low and high resolution points'
  else if (atoms%astruct%geocode =='S') then
     write(22,'(a,2x,3(1x,1pe24.17))')'surface',atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),atoms%astruct%cell_dim(3)
  else if (atoms%astruct%geocode =='P') then
     write(22,'(a,2x,3(1x,1pe24.17))')'periodic',atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),&
          atoms%astruct%cell_dim(3)
  end if
  do iat=1,atoms%astruct%nat
     write(22,'(a6,2x,3(1x,e12.5),3x)') &
          &   trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
  enddo
  do i3=0,n3  
     do i2=0,n2  
        do i1=0,n1
           if (logrid_c(i1,i2,i3))&
                &   write(22,'(a4,2x,3(1x,e10.3))') &
                &   '  g ',real(i1,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
        enddo
     enddo
  end do
  if (present(logrid_f)) then
     do i3=0,n3 
        do i2=0,n2 
           do i1=0,n1
              if (logrid_f(i1,i2,i3))&
                   &   write(22,'(a4,2x,3(1x,e10.3))') &
                   &   '  G ',real(i1,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
           enddo
        enddo
     enddo
  end if
  close(22)
END SUBROUTINE export_grids


!> Set up an array logrid(i1,i2,i3) that specifies whether the grid point
!! i1,i2,i3 is the center of a scaling function/wavelet
subroutine fill_logrid(geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,  &
      &   ntypes,iatype,rxyz,radii,rmult,hx,hy,hz,logrid)
   use module_base
   use sparsematrix_init, only: distribute_on_tasks
   implicit none
   !Arguments
   character(len=*), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
   integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nbuf,nat,ntypes
   real(gp), intent(in) :: rmult,hx,hy,hz
   integer, dimension(nat), intent(in) :: iatype
   real(gp), dimension(ntypes), intent(in) :: radii
   real(gp), dimension(3,nat), intent(in) :: rxyz
   logical, dimension(0:n1,0:n2,0:n3), intent(out) :: logrid
   !local variables
   real(kind=8), parameter :: eps_mach=1.d-12
   integer :: i1,i2,i3,iat,ml1,ml2,ml3,mu1,mu2,mu3,j1,j2,j3,i1s,i1e,i2s,i2e,i3s,i3e
   integer :: natp, isat, iiat
   real(gp) :: dx,dy2,dz2,rad,dy2pdz2,radsq
   logical :: parallel

   call f_routine(id='fill_logrid')

   !some checks
   if (geocode(1:1) /= 'F') then
      !the nbuf value makes sense only in the case of free BC
      if (nbuf /=0) then
         write(*,'(1x,a)')'ERROR: a nonzero value of nbuf is allowed only for Free BC (tails)'
         stop
      end if
   else
      !The grid spacings must be the same
      if (hx /= hy .or. hy /= hz .or. hx /= hz) then
         !        write(*,'(1x,a)')'ERROR: For Free BC the grid spacings must be the same'
      end if
   end if

   if (geocode(1:1) == 'F') then
      !$omp parallel default(none) &
      !$omp shared(nl3, nu3, nl2, nu2, nl1, nu1, logrid) &
      !$omp private(i3, i2, i1)
      !$omp do schedule(static)
      do i3=nl3,nu3 
         do i2=nl2,nu2 
            do i1=nl1,nu1
               logrid(i1,i2,i3)=.false.
            enddo
         enddo
      enddo
      !$omp end do
      !$omp end parallel
   else !
      !Special case if no atoms (homogeneous electron gas): all points are used (TD)
      if (nat == 0) then
         !$omp parallel default(none) &
         !$omp shared(n3, n2, n1, logrid) &
         !$omp private(i3, i2, i1)
         !$omp do schedule(static)
         do i3=0,n3 
            do i2=0,n2 
               do i1=0,n1
                  logrid(i1,i2,i3)=.true.
               enddo
            enddo
         enddo
         !$omp end do
         !$omp end parallel
      else
         !$omp parallel default(none) &
         !$omp shared(n3, n2, n1, logrid) &
         !$omp private(i3, i2, i1)
         !$omp do schedule(static)
         do i3=0,n3 
            do i2=0,n2 
               do i1=0,n1
                  logrid(i1,i2,i3)=.false.
               enddo
            enddo
         enddo
         !$omp end do
         !$omp end parallel
      end if
   end if

   ! MPI parallelization over the atoms, ony if there are many atoms.
   ! Maybe 200 is too low, but in this way there is a test for this feature.
   if (nat>2000) then
       call distribute_on_tasks(nat, bigdft_mpi%iproc, bigdft_mpi%nproc, natp, isat)
       parallel = .true.
   else
       natp = nat
       isat = 0
       parallel = .false.
   end if

   do iat=1,natp
      iiat = iat + isat
      rad=radii(iatype(iiat))*rmult+real(nbuf,gp)*hx
      if (rad /= 0.0_gp) then
         ml1=ceiling((rxyz(1,iiat)-rad)/hx - eps_mach)  
         ml2=ceiling((rxyz(2,iiat)-rad)/hy - eps_mach)   
         ml3=ceiling((rxyz(3,iiat)-rad)/hz - eps_mach)   
         mu1=floor((rxyz(1,iiat)+rad)/hx + eps_mach)
         mu2=floor((rxyz(2,iiat)+rad)/hy + eps_mach)
         mu3=floor((rxyz(3,iiat)+rad)/hz + eps_mach)

         !for Free BC, there must be no incoherences with the previously calculated delimiters
         if (geocode(1:1) == 'F') then
           if (ml1 < nl1) then
               write(*,'(a,i0,3x,i0)')  'ERROR: ml1 < nl1  ', ml1, nl1
               stop
           end if
           if (ml2 < nl2) then
               write(*,'(a,i0,3x,i0)')  'ERROR: ml2 < nl2  ', ml2, nl2
               stop
           end if
           if (ml3 < nl3) then
               write(*,'(a,i0,3x,i0)')  'ERROR: ml3 < nl3  ', ml3, nl3
               stop
           end if

           if (mu1 > nu1) then
               write(*,'(a,i0,3x,i0)')  'ERROR: mu1 > nu1  ', mu1, nu1
               stop
           end if
           if (mu2 > nu2) then
               write(*,'(a,i0,3x,i0)')  'ERROR: mu2 > nu2  ', mu2, nu2
               stop
           end if
           if (mu3 > nu3) then
               write(*,'(a,i0,3x,i0)')  'ERROR: mu3 > nu3  ', mu3, nu3
               stop
           end if
         end if
         i3s=max(ml3,-n3/2-1)
         i3e=min(mu3,n3+n3/2+1)
         i2s=max(ml2,-n2/2-1)
         i2e=min(mu2,n2+n2/2+1)
         i1s=max(ml1,-n1/2-1)
         i1e=min(mu1,n1+n1/2+1)
         radsq=rad**2
         !what follows works always provided the check before
         !$omp parallel default(shared) private(i3,dz2,j3,i2,dy2,j2,i1,j1,dx,dy2pdz2)
         !$omp do schedule(static,1)
         do i3=i3s,i3e
            dz2=(real(i3,gp)*hz-rxyz(3,iiat))**2-eps_mach
            if (dz2>radsq) cycle
            j3=modulo(i3,n3+1)
            do i2=i2s,i2e
               dy2=(real(i2,gp)*hy-rxyz(2,iiat))**2
               dy2pdz2=dy2+dz2
               if (dy2pdz2>radsq) cycle
               j2=modulo(i2,n2+1)
               do i1=i1s,i1e
                  j1=modulo(i1,n1+1)
                  dx=real(i1,gp)*hx-rxyz(1,iiat)
                  if (dx**2+dy2pdz2 <= radsq) then 
                     logrid(j1,j2,j3)=.true.
                  endif
               enddo
            enddo
         enddo
         !$omp enddo
         !$omp end parallel
      end if
   enddo

   if (parallel) then
       call mpiallred(logrid, mpi_lor, comm=bigdft_mpi%mpi_comm)
   end if

   call f_release_routine()

END SUBROUTINE fill_logrid


