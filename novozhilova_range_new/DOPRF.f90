module DOPRF
   use, intrinsic :: iso_c_binding

   integer(c_int) neqf, nrdf, lworkf, liworkf, aiparf(1), iparf, ioutf, ididf, itolf
   real(c_double) rtolf, atolf, rparf, ftol, artolf(1), aatolf(1), arparf(1)

   real(c_double), allocatable, target :: workf(:), yf(:)
   integer(c_int), allocatable, target :: iworkf(:)

   private allocate_arrays, deallocate_arrays
contains
   subroutine DOPRF_init(nf, ftoll)
      implicit none

      integer(c_int), intent(in) :: nf
      real(c_double), intent(in) :: ftoll

      neqf = nf
      nrdf = nf
      lworkf = 11*neqf + 8*nrdf + 21
      liworkf = nrdf + 21
      
      ftol = ftoll
      
      call allocate_arrays()
   end subroutine DOPRF_init

   subroutine allocate_arrays()
      implicit none

      integer(c_int) err_alloc

      allocate (workf(lworkf), iworkf(liworkf), yf(neqf), stat=err_alloc)

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

      deallocate (workf, iworkf, yf, stat=err_dealloc)

      if (err_dealloc /= 0) then
         print *, "deallocation error"
         pause
         stop
      end if
   end subroutine deallocate_arrays

end module DOPRF
