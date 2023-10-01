module D02PVF_D02PCF_D02PDF
   use, intrinsic :: iso_c_binding

   integer(c_int) neqp, lenwrk, method, ifailp
   real(c_double) zstart, hstart
   logical(c_bool) errass

   real(c_double), allocatable, target :: thres(:), workp(:), pgot(:), ppgot(:), pmax(:)

   private allocate_arrays, deallocate_arrays, neqp

contains

   subroutine D02PVF_D02PCF_D02PDF_init(nee)
      implicit none

      integer(c_int) nee, l

      neqp = 4*nee
      lenwrk = 32*neqp
      zstart = 0.0d0
      errass = .false.
      hstart = 0.0d0
      method = 2
      ifailp = 1

      call allocate_arrays()

      do l = 1, neqp
         thres(l) = 1.0d-8
      end do
   end subroutine D02PVF_D02PCF_D02PDF_init

   subroutine allocate_arrays()
      implicit none

      integer(c_int) err_alloc

      allocate (thres(neqp), workp(lenwrk), pgot(neqp), ppgot(neqp), pmax(neqp), stat=err_alloc)

      if (err_alloc /= 0) then
         print *, "allocation error"
         pause
         stop
      end if
   end subroutine allocate_arrays

   subroutine deallocate_arrays()
      implicit none

      integer(c_int) err_dealloc

      deallocate (thres, workp, pgot, ppgot, pmax, stat=err_dealloc)

      if (err_dealloc /= 0) then
         print *, "deallocation error"
         pause
         stop
      end if
   end subroutine deallocate_arrays

end module D02PVF_D02PCF_D02PDF
