module DOPRP
   use, intrinsic :: iso_c_binding

   integer(c_int) neqp, nrdp, lworkp, liworkp, aiparp(1), iparp, ioutp, ididp, itolp
   real(c_double) rtolp, atolp, rparp, ptol, artolp(1), aatolp(1), arparp(1)

   real(c_double), allocatable, target :: workp(:), yp(:), pex(:)
   integer(c_int), allocatable, target :: iworkp(:)

   private allocate_arrays, deallocate_arrays

contains
   subroutine DOPRP_init(ne, ptoll)
      implicit none

      integer(c_int), intent(in) :: ne
      real(c_double), intent(in) :: ptoll

      neqp = 4*ne
      nrdp = 4*ne
      lworkp = 8*neqp + 5*nrdp + 21
      liworkp = nrdp + 21

      ptol = ptoll

      call allocate_arrays()
   end subroutine DOPRP_init

   subroutine allocate_arrays()
      implicit none

      integer(c_int) err_alloc

      allocate (workp(lworkp), iworkp(liworkp), yp(neqp), pex(neqp), stat=err_alloc)

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

      deallocate (workp, iworkp, yp, pex, stat=err_dealloc)

      if (err_dealloc /= 0) then
         print *, "deallocation error"
         pause
         stop
      end if
   end subroutine deallocate_arrays

end module DOPRP
