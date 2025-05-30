!> @file
!! Include fortran file for allocation template
!! 
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if (f_nan_pad_size > 0) call togglepadding(int(0,f_long))
  if (ierror/=0) then
     !$ if(not_omp) then
     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
     !$ end if
     call f_err_throw('array ' // trim(m%array_id) // &
          & '(' // trim(yaml_toa(product(m%shape))) // &
          & '), error code '//trim(yaml_toa(ierror)),ERR_ALLOCATE)
     return
  end if
  if (bigdebug .and. any(m%shape < 0)) then
    call f_err_throw('array has suspect shape (size < 0 in at least one dimension) ' // trim(m%array_id) // &
    & '(' // trim(yaml_toa(product(m%shape))) // &
    & ')',ERR_ALLOCATE)
  end if
  if (size(shape(array))/=m%rank) then
     !$ if(not_omp) then
     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
     !$ end if
     call f_err_throw('Rank specified by f_malloc ('+yaml_toa(m%rank)// ',routine_id=' // trim(m%routine_id) // &
          & ',id=' // trim(m%array_id) // ') is not coherent with the one of the array ('+yaml_toa(size(shape(array)))//')',&
          & ERR_INVALID_MALLOC)
     return
  end if
  call pad_array(array,m%put_to_zero,m%shape,padding)
  !also fill the array with the values of the source if the address is identified in the source
