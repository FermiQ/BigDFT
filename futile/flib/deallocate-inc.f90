!> @file
!! Include fortran file for deallocation template
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  !$ if (not_omp) then
  call f_timer_interrupt(TCAT_ARRAY_ALLOCATIONS)
  !$ end if

  !here the size should be corrected with ndebug (or maybe not)
  ilsize=product(int(shape(array),kind=8))
  !retrieve the address of the first element if the size is not zero
  !iadd=int(0,kind=8)
  !if (ilsize /= int(0,kind=8)) 
  iadd=loc_arr(array)!call getlongaddress(array,iadd)

  call f_purge_database(ilsize,kind(array),iadd,info=info)

  if (associated(info)) then
     val=' '
     val=info .get. 'Type'
     select case(trim(val))
     case('SHARED')
        call bindfree(iadd)
        ierror=0
     case default
        !fallback to traditional deallocation
        deallocate(array,stat=ierror) !temporary
     end select
     call dict_free(info)
  else
     !fortran deallocation (here we should modify the calls if the array has been allocated by c)
     deallocate(array,stat=ierror)
  end if

  if (ierror/=0) then
     !$ if (not_omp) then
     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
     !$ end if
     call f_err_throw('Deallocation problem, error code '//trim(yaml_toa(ierror)),&
          ERR_DEALLOCATE)
     return
  end if

  !$ if (not_omp) then
  call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
  !$ end if
  
  
