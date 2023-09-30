module D02NVF_D02M_N
   use, intrinsic :: iso_c_binding

   integer(c_int) neqf, lwork, liwork, nrd, ifailf, it
   integer(c_int), parameter :: neq = 6, neqmax = neq, nrw = 50 + 4*neqmax, ninf = 23, nwkjac = neqmax*(neqmax + 1), &
                                maxord = 5, ny2dim = maxord + 1, maxstp = 0, mxhnil = 5
   real(c_double), parameter :: h0 = 0.0d0, hmax = 10.0d0, hmin = 1.0d-10
   integer(c_int) itask, itol, itrace, inform(ninf)
   real(c_double) const(6), rwork(nrw), wkjac(nwkjac), y(neqmax), ydot(neqmax), &
      ysave(neqmax, ny2dim), tout, atol(neqmax), rtol(neqmax), tcrit, t, xout
   logical(c_bool), parameter ::    petzld = .false.

   real(c_double), allocatable, target :: w_monitr(:, :), work(:)

   common xout, it

   private allocate_arrays, deallocate_arrays

contains

   subroutine D02NVF_D02M_N_init(ne, nt, dt, tend, ftol, f)
      implicit none

      integer(c_int), intent(in) :: ne, nt
      real(c_double), intent(in) :: tend, ftol, f(:), dt

      integer(c_int) i, it
      real(c_double) xout

      common xout, it

      neqf = 6
      nrd = 6
      lwork = 8*neqf + 5*nrd + 21
      liwork = nrd + 21

      call allocate_arrays(nt)

      !t = 0.0d0
      !tout = tend
      !itask = 1
      !itol = 1
      !rtol(1) = ftol
      !atol(1) = 0.001*ftol
      !do i = 1, 6
      !   const(i) = 0.0d0
      !end do
      !tcrit = tout
      !y(:) = f
      !itrace = -1
      !xout = dt
      !it = 2
      !ifailf = 1
   end subroutine D02NVF_D02M_N_init

   subroutine allocate_arrays(nt)
      implicit none

      integer(c_int) err_alloc, nt

      allocate (w_monitr(3, nt), work(lwork), stat=err_alloc)

      if (err_alloc /= 0) then
         print *, "allocation error"
         pause
         stop
      end if
   end subroutine allocate_arrays

   subroutine deallocate_arrays()
      implicit none

      integer(c_int) err_dealloc

      deallocate (w_monitr, work, stat=err_dealloc)

      if (err_dealloc /= 0) then
         print *, "deallocation error"
         pause
         stop
      end if
   end subroutine deallocate_arrays

end module D02NVF_D02M_N
