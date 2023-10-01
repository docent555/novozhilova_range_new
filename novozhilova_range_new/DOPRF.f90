module DOPRF
   use, intrinsic :: iso_c_binding

   integer(c_int) neqf, nrdf, lworkf, liworkf

   real(c_double), allocatable, target :: workf(:)

   private allocate_arrays, deallocate_arrays
contains
   subroutine DOPRP_init(nf)
      implicit none

      integer(c_int), intent(in) :: nf

      neqf = nf
      nrdf = nf
      lworkf = 11*neqf + 8*nrdf + 21
      liworkf = nrdf + 21
   end subroutine DOPRP_init

   subroutine allocate_arrays()
      implicit none

      integer(c_int) err_alloc

      allocate (workf(lworkf), stat=err_alloc)

      if (err_alloc /= 0) then
         print *, "allocation error"
         pause
         stop
      end if
   end subroutine allocate_arrays

   subroutine deallocate_arrays()
      use, intrinsic :: iso_c_binding
      implicit none

      integer(c_int) err_dealloc

      deallocate (workf, stat=err_dealloc)

      if (err_dealloc /= 0) then
         print *, "deallocation error"
         pause
         stop
      end if
   end subroutine deallocate_arrays

end module DOPRF
