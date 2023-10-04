module fun
   use, intrinsic :: iso_c_binding
   use ifcore
   use ifport

   integer(c_int) ne, neqp, neqf, nt, nz, freq_out, l, inharm, ndtr, METH
   real(c_double) zex_w, zex, dz, tend, dtr(2), q(3), icu(2), th(2), a(2), dcir(2), r(2), &
      f0(6), dt, pitch, f10, f20, f30, p10, p20, p30, ftol, ptol, nharm, &
      gamma, ukv, betta, betta2, betta_z, betta_z2, betta_perp, betta_perp2, gmp, w_op_w, w_op, c, e, m, b(2), w_n(2), ia(2), norm, nfac, &
      dtrb(2), dtrh_w, dtr0
   complex(c_double_complex) fp(2)
   logical(c_bool) wc, fok, lensm, btod, iatoi, inher, SQR
   character(len=6) path

   integer(c_int) breaknum(3)
   real(c_double) phitmp0(3), phitmp1(3)
   complex(c_double_complex) fc, fcomp(3)

   integer(c_int), allocatable, target :: idxre(:, :), idxim(:, :), idxp(:, :)
   complex(c_double_complex), allocatable, target :: u(:), mean(:)
   real(c_double), allocatable, target :: tax(:), zax(:), eta(:, :), etag(:, :), w(:, :), f(:, :), p(:, :), &
                                          phi(:, :), phios(:, :), wos(:, :), &
                                          cl1(:), lhs1(:), rhs1(:), cl2(:), lhs2(:), rhs2(:)

   include 'z.inc'
   include 're.inc'
   include 'im.inc'

   include 'za_22-09-23.inc'
   include 'rea_22-09-23.inc'
   include 'ima_22-09-23.inc'

   complex(c_double_complex), parameter :: ic = (0.0d0, 1.0d0)
   real(c_double), parameter :: pi = 2.0d0*dacos(0.0d0)

   private freq_out, zex_w, tend, q, i, th, a, dcir, r, f0, pitch

contains
   subroutine init()
      use D02PVF_D02PCF_D02PDF, only: D02PVF_D02PCF_D02PDF_init
      use D02NVF_D02M_N, only: D02NVF_D02M_N_init
      use DOPRF, only: DOPRF_init
      use DOPRP, only: DOPRP_init

      implicit none

      integer(c_int) ii

      call read_param()

      neqp = 4*ne
      neqf = 6

      call svch_params()

      if (lensm .eq. .true.) zex_w = betta_perp2/2.0d0/betta_z*w_op_w*zex/nharm/c
      print *, 'Zex = ', zex_w

      if (btod .eq. .true.) then
         w_n(1) = e*B(1)/(m*c)*10000.0d0
         dtr(1) = (2.0/betta_perp2)*(1.0 - (2.0*w_n(1))/(gamma*w_op_w))
         w_n(2) = e*B(2)/(m*c)*10000.0d0
         dtr(2) = (2.0/betta_perp2)*(1.0 - (2.0*w_n(2))/(gamma*w_op_w))
      end if

      if (inharm .eq. 2) then
         call norma(norm)
      else
         call norma22(norm)
      end if

      if (iatoi .eq. .true.) then
         nfac = factorial(inharm)
         icu(1) = 2.35/10000*IA(1)*(nharm**(inharm + 1)/2.0**(inharm - 1)/nfac)**2 &
                  *(Q(1)*gmp*betta**(2*inharm - 4)/gamma/betta_z)/norm
         icu(2) = 2.35/10000*IA(2)*(nharm**(inharm + 1)/2.0**(inharm - 1)/nfac)**2 &
                  *(Q(2)*gmp*betta**(2*inharm - 4)/gamma/betta_z)/norm
      end if

      print *, 'icu1 = ', icu(1), 'icu2 = ', icu(2)

      nt = tend/dt + 1
      if (nz > 0) then
         dz = zex_w/(nz - 1)
         print *, 'nz > 0'
         print *, 'nz = ', nz
         print *, 'dz = ', dz
      else
         nz = zex_w/dz + 1
         print *, 'nz == 0'
         print *, 'nz = ', nz
         print *, 'dz = ', dz
      end if

      if (METH .eq. 1) then

         call DOPRP_init(ne, ptol)
         call DOPRF_init(ne, ftol)

      else if (METH .eq. 2) then

         call D02PVF_D02PCF_D02PDF_init(ne)

         f0(1) = f10
         f0(2) = p10
         f0(3) = f20
         f0(4) = p20
         f0(5) = f30
         f0(6) = p30

         call D02NVF_D02M_N_init(6, nt, dt, tend, ftol, f0)

      else if (METH .eq. 3) then

         !call D02PVF_D02PCF_D02PDF_init(ne)
         call DOPRP_init(ne, ptol)

      end if

      call allocate_arrays()

      f(1, 1) = f10
      f(2, 1) = p10
      f(3, 1) = f20
      f(4, 1) = p20
      f(5, 1) = f30
      f(6, 1) = p30

      do ii = 1, nt
         tax(ii) = (ii - 1)*dt
      end do

      do ii = 1, nz
         zax(ii) = (ii - 1)*dz
      end do

      call calc_u(u, zex_w, nz, zax)

      do ii = 1, 2
         idxre(ii, :) = (/2*(ii - 1)*ne + 1:(2*ii - 1)*ne/)
         idxim(ii, :) = (/(2*ii - 1)*ne + 1:2*ii*ne/)
      end do

      idxp(1, :) = (/1:ne/)
      idxp(2, :) = (/ne + 1:2*ne/)

      do i = 1, ne
         p(i, 1) = dreal(cdexp(ic*(i - 1)/dble(ne)*2*pi))
         p(ne + i, 1) = dimag(cdexp(ic*(i - 1)/dble(ne)*2*pi))
         p(2*ne + i, 1) = dreal(cdexp(ic*(i - 1)/dble(ne)*2*pi))
         p(3*ne + i, 1) = dimag(cdexp(ic*(i - 1)/dble(ne)*2*pi))
      end do

      if (dtrh_w .ne. 0) then
         ndtr = (dtrb(2) - dtrb(1))/dtrh_w + 1
      else
         ndtr = 1
      end if

      dtr0 = dtrb(1)

      w(:, :) = 0

   end subroutine init

   subroutine svch_params()
      implicit none

      nharm = dble(inharm)
      gamma = 1.0 + ukv/511.0
      betta = dsqrt(1.0d0 - 1.0d0/(gamma*gamma))
      betta2 = 1.0d0 - 1.0d0/(gamma*gamma)
      betta_z = betta/dsqrt(pitch*pitch + 1.0d0)
      betta_z2 = betta2/(pitch*pitch + 1.0d0)
      betta_perp2 = betta2 - betta_z2
      w_op_w = 2*pi*w_op*1e9
      c = 29979245800.0d0
      e = 4.803e-10
      m = 9.1093837015e-28

      if (inharm .eq. 1) then
         gmp = 0.001120993533677
      else
         gmp = 0.048715056967419
      end if

   end subroutine svch_params

   subroutine init_local(path)
      implicit none

      integer(c_int) ii
      character(*) path

      call read_param_local(path)

      call svch_params()

      !if (lensm .eq. .true.) zex_w = betta_perp2/2.0d0/betta_z*w_op_w*zex/nharm/c
      !print *, 'Zex = ', zex_w

      if (btod .eq. .true.) then
         w_n(1) = e*B(1)/(m*c)*10000.0d0
         dtr(1) = (2.0/betta_perp2)*(1.0 - (2.0*w_n(1))/(gamma*w_op_w))
         w_n(2) = e*B(2)/(m*c)*10000.0d0
         dtr(2) = (2.0/betta_perp2)*(1.0 - (2.0*w_n(2))/(gamma*w_op_w))
      end if

      if (iatoi .eq. .true.) then
         nfac = factorial(inharm)
         icu(1) = 2.35/10000*IA(1)*(nharm**(inharm + 1)/2.0**(inharm - 1)/nfac)**2 &
                  *(Q(1)*gmp*betta**(2*inharm - 4)/gamma/betta_z)/norm
         icu(2) = 2.35/10000*IA(2)*(nharm**(inharm + 1)/2.0**(inharm - 1)/nfac)**2 &
                  *(Q(2)*gmp*betta**(2*inharm - 4)/gamma/betta_z)/norm
      end if

      print *, 'dtr1 = ', dtr(1), 'dtr2 = ', dtr(2)
      !print *, 'icu1 = ', icu(1), 'icu2 = ', icu(2)

      nt = tend/dt + 1
      if (nz > 0) then
         dz = zex_w/(nz - 1)
         print *, 'nz > 0'
         print *, 'nz = ', nz
         print *, 'dz = ', dz
      else
         nz = zex_w/dz + 1
         print *, 'nz == 0'
         print *, 'nz = ', nz
         print *, 'dz = ', dz
      end if

      f(1, 1) = f10
      f(2, 1) = p10
      f(3, 1) = f20
      f(4, 1) = p20
      f(5, 1) = f30
      f(6, 1) = p30

      do ii = 1, nt
         tax(ii) = (ii - 1)*dt
      end do

      do ii = 1, nz
         zax(ii) = (ii - 1)*dz
      end do

      call calc_u(u, zex_w, nz, zax)

      do ii = 1, 2
         idxre(ii, :) = (/2*(ii - 1)*ne + 1:(2*ii - 1)*ne/)
         idxim(ii, :) = (/(2*ii - 1)*ne + 1:2*ii*ne/)
      end do

      idxp(1, :) = (/1:ne/)
      idxp(2, :) = (/ne + 1:2*ne/)

      do i = 1, ne
         p(i, 1) = dreal(cdexp(ic*(i - 1)/dble(ne)*2*pi))
         p(ne + i, 1) = dimag(cdexp(ic*(i - 1)/dble(ne)*2*pi))
         p(2*ne + i, 1) = dreal(cdexp(ic*(i - 1)/dble(ne)*2*pi))
         p(3*ne + i, 1) = dimag(cdexp(ic*(i - 1)/dble(ne)*2*pi))
      end do

   end subroutine init_local

   subroutine norma(norm)
      use, intrinsic :: iso_c_binding, only: c_double, c_double_complex, c_int
      import, only:rea, ima, betta_perp2, betta_z, inharm, c, w_op_w!, omega
      implicit none

      real(c_double) :: dzz = 0.0280211, norm !v santimetrah
      complex(c_double_complex) :: u(663)

      dzz = betta_perp2/2.0d0/betta_z*w_op_w*dzz/inharm/c
      u = dcmplx(rea, ima)

      norm = sum(cdabs(u(:))*cdabs(u(:)))*dzz

      print *, 'N = ', norm

   end subroutine norma

   subroutine norma22(norm)
      use, intrinsic :: iso_c_binding, only: c_double, c_double_complex, c_int
      import, only:rea22, ima22, betta_perp2, betta_z, inharm, c, w_op_w, za22!, omega
      implicit none

      real(c_double) dzz, norm
      complex(c_double_complex) u(1123)
      integer ns

      ns = size(za22)
      dzz = sum(za22(2:1123) - za22(1:1122))/(size(za22) - 1)*0.1
      dzz = betta_perp2/2.0d0/betta_z*w_op_w*dzz/inharm/c
      u = dcmplx(rea22, ima22)

      norm = sum(cdabs(u(:))*cdabs(u(:)))*dzz

      print *, 'N = ', norm

   end subroutine norma22

   recursive function factorial(p) result(l)
      import, none
      implicit none
      integer, intent(in) :: p
      integer l
      if (p == 1) then
         l = 1
      else
         l = p*factorial(p - 1)
      end if
   end function

   function squval(zz)
      use, intrinsic :: iso_c_binding, only: c_double, c_double_complex, c_int
      import, only:zex_w, za, rea, ima
      implicit none

      real(c_double), intent(in) :: zz

      complex(c_double_complex) squval
      real(c_double) z, re, im, dz, z1
      integer(c_int) l

      dz = 0.280211
      z = zz/zex_w*185.5d0
      l = z/dz

      if (l .eq. 0) then
         l = 2
      elseif (l .ge. 662) then
         l = 662
      else
         if ((z - za(l)) .gt. 0.5*dz) l = l + 1
      end if

      z1 = za(l)
      z = z - 8.5d0

      re = rea(l - 1) + ((-rea(l - 1) + rea(l))*(dz + z - z1))/dz + ((rea(l - 1)/2.0d0 - rea(l) + &
                                                                      rea(l + 1)/2.0d0)*(z - z1)*(dz + z - z1))/dz/dz
      im = ima(l - 1) + ((-ima(l - 1) + ima(l))*(dz + z - z1))/dz + ((ima(l - 1)/2.0d0 - ima(l) + &
                                                                      ima(l + 1)/2.0d0)*(z - z1)*(dz + z - z1))/dz/dz
      !!!NE RABOTAET
      !re = ((rea(l - 1) - 2*rea(l) + rea(l + 1))*z**2)/(2.*dz**2) &
      !     + (z*(dz*(-rea(l - 1) + rea(l + 1)) - 2*(rea(l - 1) - 2*rea(l) + rea(l + 1))*z1))/(2.*dz**2) + &
      !     -(2*dz**2*rea(l) + dz*(rea(l - 1) - rea(l + 1))*z1 + (rea(l - 1) - 2*rea(l) + rea(l + 1))*z1**2)/(2.*dz**2)
      !im = ((ima(l - 1) - 2*ima(l) + ima(l + 1))*z**2)/(2.*dz**2) + &
      !     (z*(dz*(-ima(l - 1) + ima(l + 1)) - 2*(ima(l - 1) - 2*ima(l) + ima(l + 1))*z1))/(2.*dz**2) + &
      !     -(2*dz**2*ima(l) + dz*(ima(l - 1) - ima(l + 1))*z1 + (ima(l - 1) - 2*ima(l) + ima(l + 1))*z1**2)/(2.*dz**2)

      squval = dcmplx(re, im)

   end function squval

   function uval(zz)

      implicit none

      real(c_double), intent(in) :: zz

      complex(c_double_complex) uval
      real(c_double) z, re, im, d
      integer(c_int) l

      z = zz/zex_w*185.5 - 8.5
      l = (z + 8.5)/0.28021 + 1
      d = z - za(l)

      !print *, z, l, d

      if (d .gt. 0.0 .and. l /= 663) then
         re = (rea(l)*za(l + 1) - rea(l + 1)*za(l))/(za(l + 1) - za(l)) + &
              (rea(l + 1) - rea(l))/(za(l + 1) - za(l))*z
         im = (ima(l)*za(l + 1) - ima(l + 1)*za(l))/(za(l + 1) - za(l)) + &
              (ima(l + 1) - ima(l))/(za(l + 1) - za(l))*z
      else if (d .lt. 0.0 .and. l /= 1) then
         re = (rea(l - 1)*za(l) - rea(l)*za(l - 1))/(za(l) - za(l - 1)) + &
              (rea(l) - rea(l - 1))/(za(l) - za(l - 1))*z
         im = (ima(l - 1)*za(l) - ima(l)*za(l - 1))/(za(l) - za(l - 1)) + &
              (ima(l) - ima(l - 1))/(za(l) - za(l - 1))*z
      else
         re = rea(l)
         im = ima(l)
      end if

      uval = dcmplx(re, im)

   end function uval

   function squval22(zz)
      use, intrinsic :: iso_c_binding, only: c_double, c_double_complex, c_int
      import, only:zex_w, za22, rea22, ima22
      implicit none

      real(c_double), intent(in) :: zz

      complex(c_double_complex) squval22
      real(c_double) z, re, im, dz, zl
      integer(c_int) l, ns

      z = zz/zex_w*52.43d0

      l = 1
      do while (za22(l) < z)
         l = l + 1
      end do

      zl = za22(l)
      if (l .eq. 1) then
         l = 2
         zl = za22(l)
         dz = (za22(l + 1) - za22(l - 1))/2.0d0
      elseif (l .ge. 1122) then
         l = 1122
         zl = za22(l)
         dz = (za22(l + 1) - za22(l - 1))/2.0d0
      else
         if ((z - zl)/(za22(l + 1) - zl) .gt. 0.5) then
            l = l + 1
            zl = za22(l)
         end if
         dz = (za22(l + 1) - za22(l - 1))/2.0d0
      end if

      !z = z - 0.0d0

      re = rea22(l - 1) + ((-rea22(l - 1) + rea22(l))*(dz + z - zl))/dz + ((rea22(l - 1)/2.0d0 - rea22(l) + &
                                                                            rea22(l + 1)/2.0d0)*(z - zl)*(dz + z - zl))/dz/dz
      im = ima22(l - 1) + ((-ima22(l - 1) + ima22(l))*(dz + z - zl))/dz + ((ima22(l - 1)/2.0d0 - ima22(l) + &
                                                                            ima22(l + 1)/2.0d0)*(z - zl)*(dz + z - zl))/dz/dz
      !!!NE RABOTAET
      !re = ((rea22(l - 1) - 2*rea22(l) + rea22(l + 1))*z**2)/(2.*dz**2) &
      !     + (z*(dz*(-rea22(l - 1) + rea22(l + 1)) - 2*(rea22(l - 1) - 2*rea22(l) + rea22(l + 1))*zl))/(2.*dz**2) + &
      !     -(2*dz**2*rea22(l) + dz*(rea22(l - 1) - rea22(l + 1))*zl + (rea22(l - 1) - 2*rea22(l) + rea22(l + 1))*zl**2)/(2.*dz**2)
      !im = ((ima22(l - 1) - 2*ima22(l) + ima22(l + 1))*z**2)/(2.*dz**2) + &
      !     (z*(dz*(-ima22(l - 1) + ima22(l + 1)) - 2*(ima22(l - 1) - 2*ima22(l) + ima22(l + 1))*zl))/(2.*dz**2) + &
      !     -(2*dz**2*ima22(l) + dz*(ima22(l - 1) - ima22(l + 1))*zl + (ima22(l - 1) - 2*ima22(l) + ima22(l + 1))*zl**2)/(2.*dz**2)

      squval22 = dcmplx(re, im)
   end function squval22

   function uval22(zz)

      implicit none

      real(c_double), intent(in) :: zz

      complex(c_double_complex) uval22
      real(c_double) z, re, im, d
      integer(c_int) l

      z = zz/zex_w*52.43

      l = 1
      do while (za22(l) < z)
         l = l + 1
      end do
      d = z - za22(l)

      !print *, z, l, d

      if (d .gt. 0.0 .and. l /= 1123) then
         re = (rea22(l)*za22(l + 1) - rea22(l + 1)*za22(l))/(za22(l + 1) - za22(l)) + &
              (rea22(l + 1) - rea22(l))/(za22(l + 1) - za22(l))*z
         im = (ima22(l)*za22(l + 1) - ima22(l + 1)*za22(l))/(za22(l + 1) - za22(l)) + &
              (ima22(l + 1) - ima22(l))/(za22(l + 1) - za22(l))*z
      else if (d .lt. 0.0 .and. l /= 1) then
         re = (rea22(l - 1)*za22(l) - rea22(l)*za22(l - 1))/(za22(l) - za22(l - 1)) + &
              (rea22(l) - rea22(l - 1))/(za22(l) - za22(l - 1))*z
         im = (ima22(l - 1)*za22(l) - ima22(l)*za22(l - 1))/(za22(l) - za22(l - 1)) + &
              (ima22(l) - ima22(l - 1))/(za22(l) - za22(l - 1))*z
      else
         re = rea22(l)
         im = ima22(l)
      end if

      uval22 = dcmplx(re, im)

   end function uval22

   subroutine allocate_arrays()
      !use, intrinsic :: iso_c_binding
      implicit none

      integer(c_int) err_alloc

      allocate (f(6, nt), p(neqp, nz), u(nz), tax(nt), zax(nz), mean(nz), eta(2, nt), etag(2, nt), w(3, nt), &
                idxre(2, ne), idxim(2, ne), wos(3, nt), phi(3, nt), phios(3, nt), idxp(2, ne), &
                cl1(nt), lhs1(nt), rhs1(nt), cl2(nt), lhs2(nt), rhs2(nt), stat=err_alloc)

      if (err_alloc /= 0) then
         print *, "allocation error"
         pause
         stop
      end if
   end subroutine allocate_arrays

   subroutine deallocate_arrays()
      !use, intrinsic :: iso_c_binding
      implicit none

      integer(c_int) err_dealloc

      deallocate (f, p, u, tax, zax, mean, eta, etag, w, &
                  idxre, idxim, wos, phi, phios, idxp, cl1, lhs1, rhs1, cl2, lhs2, rhs2, stat=err_dealloc)

      if (err_dealloc /= 0) then
         print *, "deallocation error"
         pause
         stop
      end if
   end subroutine deallocate_arrays

   subroutine read_param_local(path) bind(c, name='read_param_local')
      use, intrinsic :: iso_c_binding
      import
      implicit none

      namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, &
         dcir1, dcir2, r1, r2, f10, f20, f30, p10, p20, p30, dt, dz, pitch, ftol, ptol, wc, fok, inharm, ukv, &
         w_op, lensm, btod, b1, b2, iatoi, ia1, ia2, &
         dtr1, dtr2, METH, SQR, nz

      real(c_double) q1, q2, q3, i1, i2, th1, th2, a1, a2, dcir1, dcir2, r1, r2, b1, b2, ia1, ia2, dtr1, dtr2
      character(*) path

      open (unit=1, file=path//'input_fortran.in', status='old', err=101)
      read (unit=1, nml=param, err=102)
      close (unit=1)

      q(1) = q1
      q(2) = q2
      q(3) = q3
      icu(1) = i1
      icu(2) = i2
      th(1) = th1
      th(2) = th2
      a(1) = a1
      a(2) = a2
      dtr(1) = dtr1
      dtr(2) = dtr2
      dcir(1) = dcir1
      dcir(2) = dcir2
      r(1) = r1
      r(2) = r2
      b(1) = b1
      b(2) = b2
      ia(1) = ia1
      ia(2) = ia2

      write (*, nml=param)

      return
101   print *, "error of file open"; pause; stop
102   print *, 'error of reading file "input_fortran.in"'; pause; stop
   end subroutine read_param_local

   subroutine read_param() bind(c, name='read_param')
      use, intrinsic :: iso_c_binding
      import
      implicit none

      namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, &
         dcir1, dcir2, r1, r2, f10, f20, f30, p10, p20, p30, dt, dz, pitch, ftol, ptol, wc, fok, inharm, ukv, &
         w_op, lensm, btod, b1, b2, iatoi, ia1, ia2, &
         dtrb, dtrh, inher, METH, SQR, nz

      real(c_double) q1, q2, q3, i1, i2, th1, th2, a1, a2, dcir1, dcir2, r1, r2, b1, b2, ia1, ia2, dtrh

      open (unit=1, file='input_fortran.in', status='old', err=101)
      read (unit=1, nml=param, err=102)
      close (unit=1)

      q(1) = q1
      q(2) = q2
      q(3) = q3
      icu(1) = i1
      icu(2) = i2
      th(1) = th1
      th(2) = th2
      a(1) = a1
      a(2) = a2
      !dtr(1) = dtr1
      !dtr(2) = dtr2
      dcir(1) = dcir1
      dcir(2) = dcir2
      r(1) = r1
      r(2) = r2
      b(1) = b1
      b(2) = b2
      ia(1) = ia1
      ia(2) = ia2
      dtrh_w = dtrh

      !write (*, nml=param)

      return
101   print *, "error of file open"; pause; stop
102   print *, 'error of reading file "input_fortran.in"'; pause; stop
   end subroutine read_param

   subroutine write_param(path) bind(c, name='write_param')
      use, intrinsic :: iso_c_binding
      import
      implicit none

      namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, &
         dcir1, dcir2, r1, r2, f10, f20, f30, p10, p20, p30, dt, dz, pitch, ftol, ptol, wc, fok, inharm, ukv, &
         w_op, lensm, btod, b1, b2, iatoi, ia1, ia2, &
         dtr1, dtr2

      character(*) path
      real(c_double) q1, q2, q3, i1, i2, th1, th2, a1, a2, dcir1, dcir2, r1, r2, b1, b2, ia1, ia2, dtr1, dtr2

      dtr1 = dtr(1)
      dtr2 = dtr(2)
      i1 = icu(1)
      i2 = icu(2)
      q1 = q(1)
      q2 = q(2)
      q3 = q(3)
      th1 = th(1)
      th2 = th(2)
      a1 = a(1)
      a2 = a(2)
      dcir1 = dcir(1)
      dcir2 = dcir(2)
      r1 = r(1)
      r2 = r(2)
      b1 = b(1)
      b2 = b(2)
      ia1 = ia(1)
      ia2 = ia(2)

      open (unit=1, file=path//'input_fortran.in', err=101)
      write (1, nml=param)
      close (unit=1)

      !write (*, nml=param)

      return
101   print *, "error of file open"; pause; stop
102   print *, 'error of reading file "input_fortran.in"'; pause; stop
   end subroutine write_param

   subroutine write_results()
      implicit none

      integer i, j

      if (wc .eq. .true.) then
         w(:, 1) = 0
         do i = 2, nt
            do j = 1, 3
               !w(j, i - 1) = dimag(log(f(2*j - 1, i)*cdexp(ic*f(2*j, i))/(f(2*j - 1, i - 1)*cdexp(ic*f(2*j, i - 1)))))/dt
               w(j, i) = (f(2*j, i) - f(2*j, i - 1))/dt
            end do
         end do
         print *, 'Frequency calculated from phase. ( WC = ', wc, ')'
      elseif (wc .eq. .false.) then
         call freq()
         print *, 'Frequency calculated from RHS. ( WC = ', wc, ')'
      end if

      phi(:, 1) = 0; 
      do i = 2, nt
         do j = 1, 3
            phi(j, i) = phi(j, i - 1) + dimag(log(f(2*j - 1, i)*cdexp(ic*f(2*j, i))/(f(2*j - 1, i - 1)*cdexp(ic*f(2*j, i - 1)))))
         end do
      end do

      breaknum(:) = 0
      fcomp(1) = f(2*1 - 1, 1)*cdexp(ic*f(2*1, 1))
      fcomp(2) = f(2*2 - 1, 1)*cdexp(ic*f(2*2, 1))
      fcomp(3) = f(2*3 - 1, 1)*cdexp(ic*f(2*3, 1))
      phitmp0(:) = datan2(dimag(fcomp(:)), dreal(fcomp(:)))
      !phitmp0(:) = datan2(dimag(f(:, 1)), dreal(f(:, 1)))
      phios(:, 1) = phitmp0(:)
      do i = 2, nt
         do j = 1, 3
            fc = f(2*j - 1, i)*cdexp(ic*f(2*j, i))
            phitmp1(j) = datan2(dimag(fc), dreal(fc))
            if ((phitmp1(j) - phitmp0(j)) .gt. pi) breaknum(j) = breaknum(j) - 1
            if ((phitmp1(j) - phitmp0(j)) .lt. -pi) breaknum(j) = breaknum(j) + 1
            phios(j, i) = phitmp1(j) + 2.*pi*breaknum(j)
            !phios(j, i) = phitmp1(j)
            phitmp0(j) = phitmp1(j)
         end do
      end do

      do i = 1, nt - 1
         do j = 1, 3
            wos(j, i) = (phios(j, i + 1) - phios(j, i))/dt
         end do
      end do

      write (*, '(/)')

      !pause

      open (3, file=path//'cl1.dat')
      do i = 1, nt
         write (3, '(5f14.6,a)') tax(i), cl1(i), lhs1(i), rhs1(i), abs(cl1(i)/lhs1(i))*100, ' %'
      end do
      close (3)

      open (3, file=path//'cl2.dat')
      do i = 1, nt
         write (3, '(5f14.6,a)') tax(i), cl2(i), lhs2(i), rhs2(i), abs(cl2(i)/lhs2(i))*100, ' %'
      end do
      close (3)

      open (1, file=path//'F.dat')
      do i = 1, nt
         !write (1, '(4e17.8)') tax(i), dabs(f(1, i)), dabs(f(3, i)), dabs(f(5, i))
         write (1, '(4f14.6)') tax(i), f(1, i), f(3, i), f(5, i)
      end do
      close (1)

      open (13, file=path//'FCMPLX.dat')
      do i = 1, nt
         fcomp(1) = f(2*1 - 1, i)*cdexp(ic*f(2*1, i))
         fcomp(2) = f(2*2 - 1, i)*cdexp(ic*f(2*2, i))
         fcomp(3) = f(2*3 - 1, i)*cdexp(ic*f(2*3, i))
         write (13, '(7f14.6)') tax(i), dreal(fcomp(1)), dimag(fcomp(1)), dreal(fcomp(2)), dimag(fcomp(2)), &
            dreal(fcomp(3)), dimag(fcomp(3))
      end do
      close (13)

      open (2, file=path//'E.dat')
      do i = 1, nt
         write (2, '(5f14.6)') tax(i), eta(1, i), etag(1, i), eta(2, i), etag(2, i)
      end do
      close (2)

      open (3, file=path//'WS.dat')
      do i = 1, nt
         write (3, '(4f14.6)') tax(i), wos(1, i), wos(2, i), wos(3, i)
      end do
      close (3)

      open (3, file=path//'W.dat')
      do i = 1, nt
         write (3, '(4f14.6)') tax(i), w(1, i), w(2, i), w(3, i)
      end do
      close (3)

      open (1, file=path//'P.dat')
      do i = 1, nt
         !write (1, '(4e17.8)') tax(i), phi(1, i), phi(2, i), phi(3, i)
         write (1, '(4f14.6)') tax(i), f(2, i), f(4, i), f(6, i)
      end do
      close (1)

      open (1, file=path//'POS.dat')
      do i = 1, nt
         write (1, '(4f14.6)') tax(i), phios(1, i), phios(2, i), phios(3, i)
      end do
      close (1)

      open (3, file=path//'WOS.dat')
      do i = 1, nt - 1
         write (3, '(4f14.6)') tax(i + 1), wos(1, i), wos(2, i), wos(3, i)
      end do
      close (3)

      !call write_param(path)

      return

101   print *, 'error of file open.'
      pause
      stop
102   print *, 'error of file reading.'
      pause
      stop
103   print *, 'error of file writing.'
      pause
      stop
   end subroutine write_results

   subroutine init_dopr_p()
      use DOPRP, only: rparp, iparp, itolp, rtolp, atolp, iworkp, workp, ptol, &
                       artolp, aatolp, arparp, aiparp
      import, only:nz, dz

      implicit none

      rparp = 0.0
      iparp = 0
      itolp = 0
      rtolp = ptol
      atolp = rtolp
      iworkp(:) = 0
      workp(:) = 0.0d0
      workp(6) = dz

      artolp(1) = rtolp
      aatolp(1) = atolp
      arparp(1) = rparp
      aiparp(1) = iparp
   end subroutine init_dopr_p

   subroutine solvef_dopr()
      use, intrinsic :: iso_c_binding
      use DOPRF, only: yf, ioutf, aatolf, arparf, itolf, workf, lworkf, liworkf, aiparf, ididf, artolf, iworkf
      implicit none

      real(c_double) rparf, t, xoutf, rtolf, atolf
      integer(c_int) iparf, itf, j

      real(c_double) pex(neqp)

      common/internf/xoutf, itf

      rparf = 0.0
      iparf = 0
      t = tax(1)
      xoutf = t
      itf = 0
      yf = f(:, 1)
      itolf = 0
      rtolf = ftol
      atolf = rtolf
      ioutf = 6

      artolf(1) = rtolf
      aatolf(1) = atolf
      arparf(1) = rparf
      aiparf(1) = iparf

      iworkf(:) = 0
      workf(:) = 0.0d0

      iworkf(5) = 6

      call dopri5_f(6, dfdt_dopr, t, yf, tend, artolf, aatolf, itolf, soloutf, ioutf, &
                    workf, lworkf, iworkf, liworkf, arparf, aiparf, ididf)

      do j = 1, neqf
         f(j, nt) = yf(j)
      end do

      call calcpex(f(:, nt), pex, cl1(nt), lhs1(nt), rhs1(nt), cl2(nt), lhs2(nt), rhs2(nt))
      eta(:, nt) = eff(pex)
      !eta(:, nt) = eff(p(:, nz))
      etag(:, nt) = pitch**2/(pitch**2 + 1)*eta(:, nt)
   end subroutine solvef_dopr

   subroutine solvep_dopr(pin, pex, c)
      use DOPRP, only: yp, iworkp, ioutp, artolp, aatolp, arparp, itolp, workp, lworkp, liworkp, aiparp, ididp

      implicit none

      real(c_double), intent(in) :: pin(:)
      real(c_double), intent(inout) :: pex(:, :)
      character(c_char), intent(in) :: c

      integer(c_int) itp
      real(c_double) z, xoutp
      common/internp/xoutp, itp

      call init_dopr_p()

      if (c .eq. 'p') then
         ioutp = 2
         z = zax(1)
         xoutp = z
         itp = 0
         yp = pin
         iworkp(5) = neqp

         call dopri5_p(neqp, dpdz_dopr, z, yp, zex_w, artolp, aatolp, itolp, soloutp, ioutp, &
                       workp, lworkp, iworkp, liworkp, arparp, aiparp, ididp)

         pex(:, nz) = yp(:)

      else
         ioutp = 0
         z = zax(1)
         xoutp = z
         itp = 0
         yp = pin
         iworkp(5) = 0

         call dopri5_p(neqp, dpdz_dopr, z, yp, zex_w, artolp, aatolp, itolp, solout_fiction, 0, &
                       workp, lworkp, iworkp, liworkp, arparp, aiparp, ididp)

         pex(:, 1) = yp(:)
      end if

   end subroutine solvep_dopr

   subroutine solvep_nag(pin, pex, c)
      use D02PVF_D02PCF_D02PDF, only: zstart, thres, method, errass, hstart, workp, lenwrk, ifailp, pgot, ppgot, pmax

      implicit none

      real(c_double), intent(in) :: pin(:)
      real(c_double), intent(inout) :: pex(:, :)
      character(c_char), intent(in) :: c

      integer(c_int) i
      real(c_double) :: zwant, zgot
      !solve eq. at t=0
      call d02pvf(neqp, zstart, p(:, 1), zex_w, ptol, thres, method, 'usual task', errass, hstart, workp, lenwrk, ifailp)

      if (c .eq. 'p') then
         do i = 1, nz - 1
            !zwant = i*dz
            zwant = zax(i + 1)
            call d02pcf(dpdz_nag, zwant, zgot, pgot, ppgot, pmax, workp, ifailp)

            if (ifailp .ne. 0) then
               write (*, *)
               write (*, 99998) 'exit d02pcf with ifail = ', ifailp, '  and z = ', zwant
               pause
               stop
            end if

            p(:, i + 1) = pgot
         end do
      else
         call d02pcf(dpdz_nag, zex_w, zgot, pex, ppgot, pmax, workp, ifailp)
      end if
99998 format(1x, a, i2, a, d12.5)
   end subroutine solvep_nag

   subroutine solvef_nag()
      use D02NVF_D02M_N, only:rtol, neqmax, ny2dim, maxord, petzld, const, tcrit, hmin, hmax, h0, maxstp, mxhnil, nwkjac, rwork, ifailf, neq, &
                               t, tout, itask, itol, atol, y, itrace, xout, it, ydot, inform, ysave, wkjac, nwkjac, w_monitr

      implicit none

      call d02nvf(neqmax, ny2dim, maxord, 'default', petzld, const, tcrit, hmin, hmax, &
                  h0, maxstp, mxhnil, 'default', rwork, ifailf)
      call d02nsf(neq, neqmax, 'a', nwkjac, rwork, ifailf)

      t = 0.0d0
      tout = tend
      itask = 1
      itol = 1
      rtol(1) = ftol
      atol(1) = 0.001*ftol
      do i = 1, 6
         const(i) = 0.0d0
      end do
      tcrit = tout
      y(:) = f(:, 1)
      itrace = -1
      xout = dt
      it = 2
      ifailf = 1

      !print *, t, tout, itask, rtol(1), atol(1), const, tcrit, y, itraice, xout, it, ifail
      !stop

      call d02nbf(neq, neqmax, t, tout, y, ydot, rwork, rtol, atol, itol, inform, dfdt_nag, ysave, &
                  ny2dim, jac, wkjac, nwkjac, monitr, itask, itrace, ifailf)

      if (ifailf .ne. 0) then
         write (*, *)
         write (*, 99998) 'exit d02nbf with ifail = ', ifailf, '  and t = ', t
         pause
         stop
      end if
99998 format(1x, a, i2, a, d12.5)
   end subroutine solvef_nag

   subroutine ode4f()
      !import, only:fp, f, ic, p, eta, etag, nz, pitch
      implicit none

      fp(1) = f(1, 1)*cdexp(ic*f(2, 1))
      fp(2) = f(3, 1)*cdexp(ic*f(4, 1))

      !solve eq. at t=0
      if (METH .eq. 1) then
         call solvep_dopr(p(:, 1), p, 'p')
      elseif (METH .eq. 2) then
         call solvep_nag(p(:, 1), p, 'p')
      elseif (METH .eq. 3) then
         call solvep_rk(p(:, 1), p, 'p')
      end if

      eta(:, 1) = eff(p(:, nz))
      etag(:, 1) = pitch**2/(pitch**2 + 1.0d0)*eta(:, 1)

      !time solver
      if (METH .eq. 1) then
         call solvef_dopr()
      elseif (METH .eq. 2) then
         call solvef_nag()
      elseif (METH .eq. 3) then
         call solvef_rk()
      end if
   end subroutine ode4f

   function eff(pex) result(eta)
      use, intrinsic :: iso_c_binding, only: c_double, c_int
      import, only:ne, idxre, idxim

      implicit none

      integer(c_int) i
      real(c_double) eta(2)
      real(c_double), intent(in) :: pex(:)

      do i = 1, 2
         eta(i) = 1 - sum(cdabs(dcmplx(pex(idxre(i, :)), pex(idxim(i, :))))**2)/ne
      end do
   end function eff

   subroutine dpdz_common(z, p, prhs)
      import :: ne, zex_w, f, ic, dtr
      implicit none

      real(c_double) z, p(*), prhs(*)

      integer(c_int) i, reidx(ne), imidx(ne)
      complex(c_double_complex) :: u
      complex(c_double_complex) s(ne), ptmp(ne)

      if (fok .eq. .false.) then
         u = dexp(-3.0d0*((z - zex_w/2)/(zex_w/2))**2)
      elseif (SQR .eq. .true.) then
         if (inharm .eq. 1) then
            u = squval22(z)
         else if (inharm .eq. 2) then
            u = squval(z)
         end if
      else
         if (inharm .eq. 1) then
            u = uval22(z)
         else if (inharm .eq. 2) then
            u = uval(z)
         end if
      end if

      do i = 1, 2
         ptmp = dcmplx(p(idxre(i, :)), p(idxim(i, :)))

         s = ic*(fp(i)*u*dconjg(ptmp)**(inharm - 1) - (dtr(i) + cdabs(ptmp)**2 - 1)*ptmp)

         prhs(idxre(i, :)) = dreal(s)
         prhs(idxim(i, :)) = dimag(s)
      end do
   end subroutine dpdz_common

   subroutine dpdz_dopr(neqp, z, p, prhs, rparp, iparp)
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex

      implicit none

      real(c_double) z, p(*), prhs(*)
      integer(c_int) neqp, iparp
      real(c_double) rparp

      call dpdz_common(z, p, prhs)
   end subroutine dpdz_dopr

   subroutine dpdz_nag(z, p, prhs)
      implicit none

      real(c_double) z, p(*), prhs(*)

      call dpdz_common(z, p, prhs)
   end subroutine dpdz_nag

   complex(c_double_complex) function xi(p, num)
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
      import, only:ne, nz, mean, u, dz, idxre, idxim, inharm

      implicit none

      integer(c_int) i, num
      real(c_double) p(:, :)

      do i = 1, nz
         mean(i) = sum(dcmplx(p(idxre(num, :), i), p(idxim(num, :), i))**inharm, 1)/ne
      end do

      mean = dconjg(u)*mean

      xi = (0.5d0*(mean(1) + mean(nz)) + sum(mean(2:(nz - 1))))*dz

   end function

   subroutine dfdt_nag(neqf, t, f, s, ires)
      implicit none

      integer(c_int) neqf, ires
      real(c_double) t, f(neqf), s(neqf)

      fp(1) = f(1)*cdexp(ic*f(2))
      fp(2) = f(3)*cdexp(ic*f(4))

      call solvep_nag(p(:, 1), p, 'p')

      call dfdt_common(neqf, t, f, s)

   end subroutine dfdt_nag

   subroutine dfdt_dopr(neqf, t, f, s, rparf, iparf)
      implicit none
      integer(c_int) :: neqf, iparf, itp
      real(c_double) t, f(neqf), s(neqf), rparf

      fp(1) = f(1)*cdexp(ic*f(2))
      fp(2) = f(3)*cdexp(ic*f(4))

      call solvep_dopr(p(:, 1), p, 'p')

      call dfdt_common(neqf, t, f, s)

   end subroutine dfdt_dopr

   subroutine dfdt_common(neqf, t, f, s)
      implicit none

      integer(c_int) :: ii, jj, iter_num = 1, time_num = 1, neqf
      real(c_double) t, f(neqf), s(neqf), &
         x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
         x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q3, &
         f1, f2, f3, phi1, phi2, phi3, a1, a2, zwant, zgot, e1, e2
      complex(c_double_complex) x1, x2
      logical bp

      x1 = xi(p, 1)
      x2 = xi(p, 2)

      x1r = dreal(x1)
      x1i = dimag(x1)
      x2r = dreal(x2)
      x2i = dimag(x2)

      f1 = f(1)
      phi1 = f(2)
      f2 = f(3)
      phi2 = f(4)
      f3 = f(5)
      phi3 = f(6)

      q31 = q(3)/q(1)
      i1 = icu(1)
      r1 = r(1)
      th1 = th(1)
      dcir1 = dcir(1)
      cos1 = dcos(phi1)
      sin1 = dsin(phi1)

      q32 = q(3)/q(2)
      i2 = icu(2)
      r2 = r(2)
      th2 = th(2)
      dcir2 = dcir(2)
      cos2 = dcos(phi2)
      sin2 = dsin(phi2)

      q3 = q(3)
      a1 = a(1)
      a2 = a(2)

      s(1) = (-nharm*f1 + i1*(-x1i*cos1 + x1r*sin1) + 2*r1*nharm*f3*dcos(phi3 - phi1 - th1))*q31
      s(2) = -2*dcir1*q3 + (i1/f1*(x1r*cos1 + x1i*sin1) + 2*r1*nharm*(f3/f1)*dsin(phi3 - phi1 - th1))*q31

      s(3) = (-nharm*f2 + i2*(-x2i*cos2 + x2r*sin2) + 2*r2*nharm*f3*dcos(phi3 - phi2 - th2))*q32
      s(4) = -2*dcir2*q3 + (i2/f2*(x2r*cos2 + x2i*sin2) + 2*r2*nharm*(f3/f2)*dsin(phi3 - phi2 - th2))*q32

      s(5) = -f3 + a1*f1*dcos(phi1 - phi3) + a2*f2*dcos(phi2 - phi3)
      s(6) = a1*f1/f3*dsin(phi1 - phi3) + a2*f2/f3*dsin(phi2 - phi3)

   end subroutine dfdt_common

   subroutine calc_u(u, zex_w, nz, zax)
      import
      implicit none

      integer(c_int), intent(in) :: nz
      real(c_double), intent(in) :: zex_w, zax(nz)
      complex(c_double_complex), intent(out) :: u(:)

      integer(c_int) i

      if (fok .eq. .false.) then
         do i = 1, nz
            u(i) = dexp(-3.0d0*((zax(i) - zex_w/2)/(zex_w/2))**2)
         end do
      elseif (SQR .eq. .true.) then
         if (inharm .eq. 1) then
            do i = 1, nz
               u(i) = squval22(zax(i))
            end do
         else if (inharm .eq. 2) then
            do i = 1, nz
               u(i) = squval(zax(i))
            end do
         end if
      else
         if (inharm .eq. 1) then
            do i = 1, nz
               u(i) = uval22(zax(i))
            end do
         else if (inharm .eq. 2) then
            do i = 1, nz
               u(i) = uval(zax(i))
            end do
         end if
      end if

      open (1, file='struc.dat')
      do i = 1, nz
         if (inharm == 2) then
            write (1, '(3f14.6)') zax(i)/zax(nz)*185.5 - 8.5, dreal(u(i)), dimag(u(i))
         else
            write (1, '(3f14.6)') zax(i)/zax(nz)*52.43, dreal(u(i)), dimag(u(i))
         end if
      end do
      close (1)
      !stop

   end subroutine

   subroutine jac(neq, t, y, h, d, pp)
      implicit none
!     .. scalar arguments ..
      double precision d, h, t
      integer neq
!     .. array arguments ..
      double precision pp(neq, neq), y(neq)
!     .. local scalars ..
      double precision hxd
!     .. executable statements ..

      real(c_double) x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
         x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q333, &
         f1, f2, f3, phi1, phi2, phi3, a1, a2, ph311, ph322, phi13, phi23, &
         sin311, sin322, cos311, cos322, cos13, sin13, cos23, sin23, phi131, phi232
      complex(c_double_complex) x1, x2

      hxd = h*d

      x1 = xi(p, 1)
      x2 = xi(p, 2)

      x1r = dreal(x1)
      x1i = dimag(x1)
      x2r = dreal(x2)
      x2i = dimag(x2)

      f1 = y(1)
      phi1 = y(2)
      f2 = y(3)
      phi2 = y(4)
      f3 = y(5)
      phi3 = y(6)

      q31 = q(3)/q(1)
      i1 = icu(1)
      r1 = r(1)
      th1 = th(1)
      dcir1 = dcir(1)
      cos1 = dcos(phi1)
      sin1 = dsin(phi1)

      q32 = q(3)/q(2)
      i2 = icu(2)
      r2 = r(2)
      th2 = th(2)
      dcir2 = dcir(2)
      cos2 = dcos(phi2)
      sin2 = dsin(phi2)

      ph311 = phi3 - phi1 - th1
      ph322 = phi3 - phi2 - th2
      phi13 = phi1 - phi3
      phi23 = phi2 - phi3
      phi131 = phi1 - phi3 + th1
      phi232 = phi2 - phi3 + th2

      cos311 = dcos(ph311)
      sin311 = dsin(ph311)
      cos322 = dcos(ph322)
      sin322 = dsin(ph322)
      cos13 = dcos(phi13)
      sin13 = dsin(phi13)
      cos23 = dcos(phi23)
      sin23 = dsin(phi23)

      !pp(1, 1) = 1.0d0 + hxd*q31*nharm
      !pp(1, 2) = -hxd*q31*(i1*(x1i*sin1 + x1r*cos1) + 2.0d0*nharm*r1*f3*sin311)
      !!pp(1,3) = 0
      !!pp(1,4) = 0
      !pp(1, 5) = -hxd*2.0d0*nharm*r1*q31*cos311
      !pp(1, 6) = hxd*2.0d0*nharm*r1*q31*f3*sin311
      !
      !pp(2, 1) = hxd*q31/(f1**2)*(i1*(x1r*cos1 + x1i*sin1) + 2.0d0*nharm*r1*f3*sin311)
      !pp(2, 2) = 1.0d0 - hxd*q31/f1*(i1*(-x1r*sin1 + x1i*cos1) - 2.0d0*nharm*r1*f3*cos311)
      !!pp(2,3) = 0
      !!pp(2,4) = 0
      !pp(2, 5) = -hxd*2.0d0*nharm*r1*q31/f1*sin311
      !pp(2, 6) = -hxd*2.0d0*nharm*r1*q31*f3/f1*cos311
      !
      !!pp(3,1) = 0
      !!pp(3,2) = 0
      !pp(3, 3) = 1.0d0 + hxd*q32*nharm
      !pp(3, 4) = -hxd*q32*(i2*(x2i*sin2 + x2r*cos2) + 2.0d0*nharm*r2*f3*sin322)
      !pp(3, 5) = -hxd*2.0d0*nharm*r2*q32*cos322
      !pp(3, 6) = hxd*2.0d0*nharm*r2*q32*f3*sin322
      !
      !!pp(4,1) = 0
      !!pp(4,2) = 0
      !pp(4, 3) = hxd*q32/(f2**2)*(i2*(x2r*cos2 + x2i*sin2) + 2.0d0*nharm*r2*f3*sin322)
      !pp(4, 4) = 1.0d0 - hxd*q32/f2*(i2*(-x2r*sin2 + x2i*cos2) - 2.0d0*nharm*r2*f3*cos322)
      !pp(4, 5) = -hxd*2.0d0*nharm*r2*q32/f2*sin322
      !pp(4, 6) = -hxd*2.0d0*nharm*r2*q32*f3/f2*cos322
      !
      !pp(5, 1) = -hxd*a1*cos13
      !pp(5, 2) = hxd*a1*f1*sin13
      !pp(5, 3) = -hxd*a2*cos23
      !pp(5, 4) = hxd*a2*f2*sin23
      !pp(5, 5) = 1.0d0 + hxd
      !pp(5, 6) = -hxd*(a1*f1*sin13 + a2*f2*sin23)
      !
      !pp(6, 1) = -hxd*a1/f3*sin13
      !pp(6, 2) = -hxd*a1*f1/f3*cos13
      !pp(6, 3) = -hxd*a2/f3*sin23
      !pp(6, 4) = -hxd*a2*f2/f3*cos23
      !pp(6, 5) = hxd/(f3**2)*(a1*f1*sin13 + a2*f2*sin23)
      !pp(6, 6) = 1.0d0 + hxd/f3*(a1*f1*cos13 + a2*f2*cos23)

      pp(1, 1) = hxd*nharm*q31 + 1.0D0
      pp(1, 2) = -hxd*q31*(i1*(x1r*cos1 + x1i*sin1) - f3*nharm*r1*sin(phi131)*2.0D0)
      !pp(1, 3) = 0.0D0
      !pp(1, 4) = 0.0D0
      pp(1, 5) = hxd*nharm*q31*r1*cos(phi131)*(-2.0D0)
      pp(1, 6) = f3*hxd*nharm*q31*r1*sin(phi131)*(-2.0D0)

      pp(2, 1) = hxd*q31*(1.0D0/f1**2*i1*(x1r*cos1 + x1i*sin1) - 1.0D0/f1**2*f3*nharm*r1*sin(phi131)*2.0D0)
      pp(2, 2) = -hxd*q31*((i1*(x1i*cos1 - x1r*sin1))/f1 - (f3*nharm*r1*cos(phi131)*2.0D0)/f1) + 1.0D0
      !pp(2, 3) = 0.0D0
      !pp(2, 4) = 0.0D0
      pp(2, 5) = (hxd*nharm*q31*r1*sin(phi131)*2.0D0)/f1
      pp(2, 6) = (f3*hxd*nharm*q31*r1*cos(phi131)*(-2.0D0))/f1

      !pp(3, 1) = 0.0D0
      !pp(3, 2) = 0.0D0
      pp(3, 3) = hxd*nharm*q32 + 1.0D0
      pp(3, 4) = -hxd*q32*(i2*(x2r*cos2 + x2i*sin2) - f3*nharm*r2*sin(phi232)*2.0D0)
      pp(3, 5) = hxd*nharm*q32*r2*cos(phi232)*(-2.0D0)
      pp(3, 6) = f3*hxd*nharm*q32*r2*sin(phi232)*(-2.0D0)

      !pp(4, 1) = 0.0D0
      !pp(4, 2) = 0.0D0
      pp(4, 3) = hxd*q32*(1.0D0/f2**2*i2*(x2r*cos2 + x2i*sin2) - 1.0D0/f2**2*f3*nharm*r2*sin(phi232)*2.0D0)
      pp(4, 4) = -hxd*q32*((i2*(x2i*cos2 - x2r*sin2))/f2 - (f3*nharm*r2*cos(phi232)*2.0D0)/f2) + 1.0D0
      pp(4, 5) = (hxd*nharm*q32*r2*sin(phi232)*2.0D0)/f2
      pp(4, 6) = (f3*hxd*nharm*q32*r2*cos(phi232)*(-2.0D0))/f2

      pp(5, 1) = -a1*hxd*cos13
      pp(5, 2) = a1*f1*hxd*sin13
      pp(5, 3) = -a2*hxd*cos23
      pp(5, 4) = a2*f2*hxd*sin23
      pp(5, 5) = hxd + 1.0D0
      pp(5, 6) = -hxd*(a1*f1*sin13 + a2*f2*sin23)

      pp(6, 1) = -(a1*hxd*sin13)/f3
      pp(6, 2) = -(a1*f1*hxd*cos13)/f3
      pp(6, 3) = -(a2*hxd*sin23)/f3
      pp(6, 4) = -(a2*f2*hxd*cos23)/f3
      pp(6, 5) = hxd*(a1*f1*1.0D0/f3**2*sin13 + a2*f2*1.0D0/f3**2*sin23)
      pp(6, 6) = hxd*((a1*f1*cos13)/f3 + (a2*f2*cos23)/f3) + 1.0D0

      return
   end

   subroutine monitr(n, nmax, t, hlast, h, y, ydot, ysave, r, acor, imon, inln, hmin, hmxi, nqu)
      use D02NVF_D02M_N, only: w_monitr
      implicit none
      !..parameters..
      integer nout, it, j
      parameter(nout=6)
      integer ny2dim
      parameter(ny2dim=6)
      !..scalar arguments..
      double precision h, hlast, hmin, hmxi, t
      integer imon, inln, n, nmax, nqu
      !..array arguments..
      double precision acor(nmax, 2), r(n), y(n), ydot(n), ysave(nmax, *)
      !..scalars in common..
      double precision xout
      !..local scalars..
      integer i, ifail
      logical(4) pressed
      character(1) key
      integer(c_int), parameter :: esc = 27

      real(c_double) pex(neqp)

      !..external subroutines..
      external d02xkf
      !..common blocks..
      common xout, it
      !..executable statements..
      if (imon .ne. 1) return
20    if (.not. (t - hlast .lt. xout .and. xout .le. t)) return
      ifail = 1
      !c1 interpolation
      call d02xkf(xout, r, n, ysave, nmax, ny2dim, acor(1, 2), n, t, nqu, hlast, h, ifail)

      if (ifail .ne. 0) then
         imon = -2
      else
         !write (nout, 99999) xout, (r(i), i=1, n)
         call calcpex(r, pex, cl1(it), lhs1(it), rhs1(it), cl2(it), lhs2(it), rhs2(it))
         !call zs(r, pex, cl1(it), lhs1(it), rhs1(it), cl2(it), lhs2(it), rhs2(it))
         eta(:, it) = eff(pex)
         etag(:, it) = pitch**2/(pitch**2 + 1)*eta(:, it)
         write (*, '(a,f8.3,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f9.5,a,f9.5,a,f9.5,a,f5.3,a,f5.3,a,\,a)') 't =', xout, &
            '  |f1| =', r(1), '  |f2| =', r(3), '  |f3| =', r(5), '  e1 =', eta(1, it), '  e2 =', eta(2, it), &
            '  w1 = ', ydot(2), '  w2 = ', ydot(4), '  w3 = ', ydot(6), '  c1 = ', dabs(cl1(it)/rhs1(it))*100, '%  c2 = ', dabs(cl2(it)/rhs2(it)*100), '%', char(13)
         do j = 1, n
            f(j, it) = r(j)
         end do
         do j = 1, 3
            w_monitr(j, it) = ydot(2*j)
            w(j, it) = w_monitr(j, it)
         end do
         xout = xout + dt
         nt = it
         it = it + 1
         pressed = peekcharqq()
         if (pressed) then
            key = getcharqq()
            if (ichar(key) .eq. esc) then
               write (*, '(/,a)') 'quit?'
               key = getcharqq()
               if (ichar(key) .eq. 121 .or. ichar(key) .eq. 89) then
                  nt = it - 1
                  call write_results()
                  !imon = -2
                  !return
               end if
            end if
         end if
         if (xout .lt. tend) go to 20
      end if

      return

99999 format(1x, f8.3, 6(f13.5, 2x))
   end subroutine monitr

   subroutine calcpex(fpex, pex, c1, lhs1, rhs1, c2, lhs2, rhs2)

      implicit none

      real(8), intent(out) :: c1, lhs1, rhs1, c2, lhs2, rhs2

      integer(c_int) i, idx(ne)
      real(c_double) :: zwant, zgot, pex(neqp, 1), fpex(6), p2ex_mean(2), p20_mean(2)
      complex(8) ptmp(ne)

      fp(1) = fpex(1)*cdexp(ic*fpex(2))
      fp(2) = fpex(3)*cdexp(ic*fpex(4))

      if (METH .eq. 1) then
         call solvep_dopr(p(:, 1), pex, 'a')
      elseif (METH .eq. 2) then
         call solvep_nag(p(:, 1), pex, 'a')
      elseif (METH .eq. 3) then
         !call solvep_rk(p(:, 1), pex, 'a')
         !call solvep_nag(p(:, 1), pex, 'a')
         call solvep_dopr(p(:, 1), pex, 'a')
      end if

      do i = 1, 2
         idx = idxp(i, :)

         ptmp(:) = dcmplx(p(idxre(i, :), 1), p(idxim(i, :), 1))
         p20_mean(i) = sum(cdabs(ptmp(:)*cdabs(ptmp(:))))/ne

         ptmp(:) = dcmplx(pex(idxre(i, :), 1), pex(idxim(i, :), 1))
         p2ex_mean(i) = sum(cdabs(ptmp(:)*cdabs(ptmp(:))))/ne
      end do

      lhs1 = 2.0*nharm*fpex(1)**2 - 4.0*nharm*r(1)*fpex(5)*fpex(1)*dcos(th(1) - fpex(6) + fpex(2))
      rhs1 = -icu(1)*(p2ex_mean(1) - p20_mean(1))
      c1 = lhs1 - rhs1

      lhs2 = 2.0*nharm*fpex(3)**2 - 4.0*nharm*r(2)*fpex(5)*fpex(3)*dcos(th(2) - fpex(6) + fpex(4))
      rhs2 = -icu(2)*(p2ex_mean(2) - p20_mean(2))
      c2 = lhs2 - rhs2

   end subroutine calcpex

   subroutine freq()

      implicit none

      integer(c_int) :: ii, iparf, aiparp(1), itp, j
      real(c_double) t, z, zwant, zgot, &
         x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
         x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q3, &
         f1, f2, f3, phi1, phi2, phi3, a1, a2
      complex(c_double_complex) x1, x2

      do ii = 2, nt

         fp(1) = f(1, ii)*exp(ic*f(2, ii))
         fp(2) = f(3, ii)*exp(ic*f(4, ii))

         if (METH .eq. 1) then
            call solvep_dopr(p(:, 1), p, 'p')
         elseif (METH .eq. 2) then
            call solvep_nag(p(:, 1), p, 'p')
         elseif (METH .eq. 3) then
            call solvep_rk(p(:, 1), p, 'p')
         end if

         !x1 = xi(p(1:2*ne, :), 1)
         !x2 = xi(p(2*ne + 1:4*ne, :), 1)
         x1 = xi(p, 1)
         x2 = xi(p, 2)

         x1r = dreal(x1)
         x1i = dimag(x1)
         x2r = dreal(x2)
         x2i = dimag(x2)

         f1 = f(1, ii)
         phi1 = f(2, ii)
         f2 = f(3, ii)
         phi2 = f(4, ii)
         f3 = f(5, ii)
         phi3 = f(6, ii)

         q31 = q(3)/q(1)
         i1 = icu(1)
         r1 = r(1)
         th1 = th(1)
         dcir1 = dcir(1)
         cos1 = dcos(phi1)
         sin1 = dsin(phi1)

         q32 = q(3)/q(2)
         i2 = icu(2)
         r2 = r(2)
         th2 = th(2)
         dcir2 = dcir(2)
         cos2 = dcos(phi2)
         sin2 = dsin(phi2)

         q3 = q(3)
         a1 = a(1)
         a2 = a(2)

         !s(1) = (-nharm*f1 + i1*(-x1i*cos1 + x1r*sin1) + 2*r1*nharm*f3*dcos(phi3 - phi1 - th1))*q31
         w(1, ii) = -2*dcir1*q3 + (i1/f1*(x1r*cos1 + x1i*sin1) + 2*r1*nharm*(f3/f1)*dsin(phi3 - phi1 - th1))*q31

         !s(3) = (-nharm*f2 + i2*(-x2i*cos2 + x2r*sin2) + 2*r2*nharm*f3*dcos(phi3 - phi2 - th2))*q32
         w(2, ii) = -2*dcir2*q3 + (i2/f2*(x2r*cos2 + x2i*sin2) + 2*r2*nharm*(f3/f2)*dsin(phi3 - phi2 - th2))*q32

         !s(5) = -f3 + a1*f1*dcos(phi1 - phi3) + a2*f2*dcos(phi2 - phi3)
         w(3, ii) = a1*f1/f3*dsin(phi1 - phi3) + a2*f2/f3*dsin(phi2 - phi3)
      end do
   end subroutine freq

   subroutine solout_fiction
   end subroutine solout_fiction

   subroutine soloutp(nr, xold, x, y, n, con, icomp, nd, rparp, iparp, irtrn)
      implicit none

      interface
         function contd5_p(ii, x, con, icomp, nd)
            implicit double precision(a - h, o - z)
            dimension con(5*nd), icomp(nd)
         end
      end interface

      integer(c_int) nr, n, nd, icomp(nd), iparp, irtrn, j, itp
      real(c_double) xold, x, con(5*nd), rparp, y(neqp), xoutp
      logical(4) pressed
      character(1) key
      integer(c_int), parameter :: esc = 27
      common/internp/xoutp, itp

      if (nr .eq. 1) then
         itp = 1
         do j = 1, neqp
            p(j, itp) = y(j)
         end do
         xoutp = x + dz
      else
10       continue
         if (x .ge. xoutp) then
            itp = itp + 1
            do j = 1, neqp
               p(j, itp) = contd5_p(j, xoutp, con, icomp, nd)
            end do
            xoutp = xoutp + dz
            goto 10
         end if
      end if
      return
   end subroutine soloutp

   subroutine soloutf(nr, xold, x, y, n, con, icomp, nd, rparf, iparf, irtrn)
      implicit none

      interface
         function contd5_f(ii, x, con, icomp, nd)
            implicit double precision(a - h, o - z)
            dimension con(5*nd), icomp(nd)
         end
      end interface

      integer(c_int) nr, n, nd, icomp(nd), iparf, irtrn, j, itf
      real(c_double) xold, x, con(5*nd), rparf, y(neqf), xoutf, pex(neqp), yy(neqf), wcur(3)
      logical(4) pressed
      character(1) key
      integer(c_int), parameter :: esc = 27
      common/internf/xoutf, itf

      if (nr .eq. 1) then
         itf = 1
         do j = 1, neqf
            f(j, itf) = y(j)
         end do
         call calcpex(y, pex, cl1(itf), lhs1(itf), rhs1(itf), cl2(itf), lhs2(itf), rhs2(itf))
         eta(:, itf) = eff(pex)
         !eta(:, itf) = eff(p(:, nz))
         etag(:, itf) = pitch**2/(pitch**2 + 1)*eta(:, itf)
         !write (*, '(a,f10.5,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f6.3,a,f6.3,\,a,a)') 'Time = ', xoutf, &
         write (*, '(a,f8.3,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f9.5,a,f9.5,a,f9.5,a,f5.3,a,f5.3,a,\,a)') 't =', xoutf, &
            '  |F1| = ', abs(f(1, itf)), '  |F2| = ', abs(f(3, itf)), &
            '  |F3| = ', abs(f(5, itf)), '  Eff1 = ', eta(1, itf), '  Eff2 = ', eta(2, itf), &
            '  w1 = ', 0, '  w2 = ', 0, '  w3 = ', 0, &
            !'  ph1 = ', f(2, itf), '  ph2 = ', f(4, itf), '  ph3 = ', f(6, itf), &
            '  c1 = ', abs(cl1(itf)/rhs1(itf))*100, ' %  c2 = ', abs(cl2(itf)/rhs2(itf))*100, ' %', char(13)
         xoutf = x + dt
      else
10       continue
         if (x .ge. xoutf) then
            itf = itf + 1
            do j = 1, neqf
               f(j, itf) = contd5_f(j, xoutf, con, icomp, nd)
            end do
            call calcpex(y, pex, cl1(itf), lhs1(itf), rhs1(itf), cl2(itf), lhs2(itf), rhs2(itf))
            eta(:, itf) = eff(pex)
            !eta(:, itf) = eff(p(:, nz))
            etag(:, itf) = pitch**2/(pitch**2 + 1)*eta(:, itf)
            do j = 1, 3
               wcur(j) = (f(2*j, itf) - f(2*j, itf - 1))/dt
            end do
            !write (*, '(a,f10.5,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f6.3,a,f6.3,\,a,a)') 'Time = ', xoutf, &
            write (*, '(a,f8.3,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f9.5,a,f9.5,a,f9.5,a,f5.3,a,f5.3,a,\,a)') 't =', xoutf, &
               '  |F1| = ', abs(f(1, itf)), '  |F2| = ', abs(f(3, itf)), &
               '  |F3| = ', abs(f(5, itf)), '  Eff1 = ', eta(1, itf), '  Eff2 = ', eta(2, itf), &
               '  w1 = ', wcur(1), '  w2 = ', wcur(2), '  w3 = ', wcur(3), &
               '  c1 = ', dabs(cl1(itf)/rhs1(itf))*100, ' %  c2 = ', dabs(cl2(itf)/rhs2(itf))*100, ' %', char(13)
            xoutf = xoutf + dt
            goto 10
         end if
      end if

      pressed = peekcharqq()
      if (pressed) then
         key = getcharqq()
         if (ichar(key) .eq. esc) then
            write (*, '(/,a)') 'Quit?'
            key = getcharqq()
            if (ichar(key) .eq. 121 .or. ichar(key) .eq. 89) then
               nt = itf
               irtrn = -1
               !return
            end if
         end if
      end if
      return
   end subroutine soloutf

   subroutine solvef_rk()
      implicit none

      !integer(c_int) nt, neq
      real(c_double) h, t0, t

      call ode4f_rk(dfdt_rk, f, neqf, nt, 0.0d0, dt)
   end subroutine solvef_rk

   function dfdt_rk(t, y) result(s)
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double) t, y(:), s(size(y))

      call dfdt_common(neqf, t, y, s)
   end function dfdt_rk

   subroutine ode4f_rk(dydt, y, neq, nt, t0, h)
      import
      implicit none

      interface
         function dydt(t, y) result(s)
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double) t, y(:), s(size(y))
         end function dydt
      end interface

      integer(c_int) nt, neq, i, j, itf
      real(c_double) h, t0, t
      real(c_double), intent(inout) :: y(:, :)
      real(c_double) s1(size(y, 1)), s2(size(y, 1)), s3(size(y, 1)), s4(size(y, 1)), v(size(y, 1))
      real(c_double) pex(neqp)
      logical(4) pressed
      character(1) key
      integer(c_int), parameter :: esc = 27

      !!solve eq. at t=0
      !fp(1) = y(1, 1)*cdexp(ic*y(2, 1))
      !fp(2) = y(3, 1)*cdexp(ic*y(4, 1))
      !
      !call ode4p(dpdz_common, p, ne, nz, 0.0d0, dz)
      !
      !eta(:, 1) = eff(p(:, nz))
      !etag(:, 1) = pitch**2/(pitch**2 + 1)*eta(:, 1)

      !open (1, file='test.dat')
      !do j = 1, nz
      !    write (1, '(5f17.8)') (j-1)*dz, p(14, j), p(ne+14, j), p(128, j), p(ne+128, j)
      !end do
      !close (1)
      !stop

      do i = 1, nt - 1
         v = y(:, i)
         t = t0 + (i - 1)*h
         s1(:) = dydt(t, v)
         s2(:) = dydt(t + h/2, v + h*s1(:)/2)
         s3(:) = dydt(t + h/2, v + h*s2(:)/2)
         s4(:) = dydt(t + h, v + h*s3(:))
         y(:, i + 1) = v + h*(s1(:) + 2.0d0*s2(:) + 2.0d0*s3(:) + s4(:))/6.0d0

         eta(:, i + 1) = eff(p(:, nz))
         etag(:, i + 1) = pitch**2/(pitch**2 + 1)*eta(:, i + 1)

         itf = i + 1
         call calcpex(y(:, i + 1), pex, cl1(itf), lhs1(itf), rhs1(itf), cl2(itf), lhs2(itf), rhs2(itf))
         eta(:, itf) = eff(pex)
         !eta(:, itf) = eff(p(:, nz))
         etag(:, itf) = pitch**2/(pitch**2 + 1)*eta(:, itf)
         w(1, itf) = (s1(2) + 2.0d0*s2(2) + 2.0d0*s3(2) + s4(2))/6.0d0
         w(2, itf) = (s1(4) + 2.0d0*s2(4) + 2.0d0*s3(4) + s4(4))/6.0d0
         w(3, itf) = (s1(6) + 2.0d0*s2(6) + 2.0d0*s3(6) + s4(6))/6.0d0

         !write (*, '(a,f12.7,a,f10.7,a,f10.7,a,f10.7,a,f10.7,a,f10.7,a,f5.3,a,f5.3,a,\,a)') 'Time = ', t + h, '   |F1| = ', abs(y(1, i + 1)), '   |F2| = ', abs(y(3, i + 1)), '   |F3| = ', abs(y(5, i + 1)), &
         !'   Eff1 = ', eta(1, i + 1), '   Eff2 = ', eta(2, i + 1), '  c1 = ', dabs(cl1(itf)/rhs1(itf))*100, ' %  c2 = ', dabs(cl2(itf)/rhs2(itf))*100, ' %', char(13)
         write (*, '(a,f8.3,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f9.5,a,f9.5,a,f9.5,a,f5.3,a,f5.3,a,\,a)') 't =', t + h, &
            '  |f1| =', abs(y(1, i + 1)), '  |f2| =', abs(y(3, i + 1)), '  |f3| =', abs(y(5, i + 1)), '  e1 =', eta(1, i + 1), '  e2 =', eta(2, i + 1), &
            '  w1 = ', w(1, itf), '  w2 = ', w(2, itf), '  w3 = ', w(3, itf), '  c1 = ', dabs(cl1(itf)/rhs1(itf))*100, '%  c2 = ', dabs(cl2(itf)/rhs2(itf)*100), '%', char(13)

         pressed = peekcharqq()
         if (pressed) then
            key = getcharqq()
            if (ichar(key) .eq. esc) then
               write (*, '(/,a)') 'Quit?'
               key = getcharqq()
               if (ichar(key) .eq. 121 .or. ichar(key) .eq. 89) then
                  nt = i + 1
                  return
                  !open (1, file='F.dat')
                  !do j = 1, i + 1
                  !    write (1, '(4e17.8)') (j - 1)*h, abs(y(1, j)), abs(y(2, j)), abs(y(3, j))
                  !end do
                  !close (2)
                  !open (2, file='E.dat')
                  !do j = 1, i + 1
                  !    write (2, '(5e17.8)') (j - 1)*h, eta(1, j), etag(1, j), eta(2, j), etag(2, j)
                  !end do
                  !close (2)
                  !stop
               end if
            end if
         end if
      end do
   end subroutine ode4f_rk

   subroutine solvep_rk(pin, pex, c)
      implicit none

      real(c_double), intent(in) :: pin(:)
      real(c_double), intent(inout) :: pex(:, :)
      character(c_char), intent(in) :: c

      pex(:, 1) = pin

      if (c .eq. 'p') then
         call ode4p_rk(dpdz_rk, pex, neqp, nz, 0.0d0, dz)
      else
         call ode4p_rk(dpdz_rk, p, neqp, nz, 0.0d0, dz)
         pex(:, 1) = p(:, nz)
      end if
   end subroutine solvep_rk

   subroutine ode4p_rk(dydt, y, neq, nz, z0, h)!, params)
      import
      implicit none

      interface
         subroutine dydt(z, p, prhs)
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double), intent(in) :: z, p(*)
            real(c_double), intent(out) :: prhs(*)
         end subroutine dydt
      end interface

      integer(c_int) nz, neq, i
      real(c_double) h, z0, z
      real(c_double), intent(inout) :: y(:, :)
      real(c_double) s1(size(y, 1)), s2(size(y, 1)), s3(size(y, 1)), s4(size(y, 1)), v(size(y, 1))

      do i = 1, nz - 1
         v = y(:, i)
         z = z0 + (i - 1)*h
         call dydt(z, v, s1)
         call dydt(z + h/2, v + h*s1(:)/2, s2)
         call dydt(z + h/2, v + h*s2(:)/2, s3)
         call dydt(z + h, v + h*s3(:), s4)
         y(:, i + 1) = v + h*(s1(:) + 2*s2(:) + 2*s3(:) + s4(:))/6
      end do
   end subroutine ode4p_rk

   subroutine dpdz_rk(z, p, prhs)
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), intent(in) :: z, p(*)
      real(c_double), intent(out) :: prhs(*)

      call dpdz_common(z, p, prhs)

   end subroutine dpdz_rk

end module fun
