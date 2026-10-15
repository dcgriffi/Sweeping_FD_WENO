PROGRAM FDWENO4EULER1D
   USE weno5
   use MPI
   IMPLICIT NONE
!********************************************************
! This is a parallel MPI code for 3rd order RK + 5th order
! finite difference WENO for the 2D Euler equation.
!                     Xiangxiong Zhang  Brown University
!                                         July 2011
!********************************************************
   INTEGER, PARAMETER:: total = 80, cut_x = 4, cut_y = 8, nx = 5*total/cut_x,& ! cut_x should be less than or equal to 5
  & ny = 5*total/cut_y, ierror = 0, mg = 0, step = 2
   !t_final = 0.6
   INTEGER, PARAMETER :: total_x = nx*cut_x, total_y = ny*cut_y

   !time_terminal = 0.6
   REAL, PARAMETER:: time_terminal = 0.6, gamma = 1.2, gamma1 = gamma - 1., eps = 1d-13, eps2 = 1d-8
   REAL, PARAMETER:: q_const = 50., t_const = 50., k_const = 2566.4

   REAL, DIMENSION(5, -2:nx + 3, -2:ny + 3):: u, u0, RHS, fp_flux, fn_flux, fp, fn, gp_flux, gn_flux, gp, gn
   REAL, DIMENSION(-2:nx + 3):: x
   REAL, DIMENSION(-2:ny + 3):: y
   REAL, DIMENSION(-2:nx + 3, -2:ny + 3)::alpha
   integer:: rktag, kcount
   real:: dx, dy, dt, time, alphamax, xlen, ylen, source_alpha
   integer:: ig, jg, i, j, k, kk, m, den_count = 0, pre_count = 0, postag

   REAL(8):: utemp(5, 1:total_x, 1:total_y)

   INTEGER:: mpirank, mpierror, mpisize, position_x, position_y, destination(8), source(8)&
  &, mpitag, mpisource, recv_status(MPI_STATUS_SIZE), send_status(MPI_STATUS_SIZE), mpides&
  &, send_request(4), recv_request(4), color_x, color_y
   REAL:: time_start, time_end
   real:: xm1, ss, dri, uri, vri, pri, ach2, ple, ule, dle, vle

   !Adding this for dynamics output storage
   character(len=256) :: output_dir, err_file_path, den_file_path
   character(len=256) :: pre_file_path, rate_file_path
   integer :: ios

   !Adding sweeping types
   integer::  x_iter, y_iter, iter
   INTEGER, PARAMETER :: use_limiter = 2, swp_num = 20000
   Real(8), PARAMETER:: nSteps_x = total_x, nsteps_y = total_y

   INTEGER, PARAMETER :: local_size = nx*ny
   REAL(8), DIMENSION(1:total_x, 1:total_y):: weights
   integer, allocatable :: sweepI_i(:), sweepI_j(:)
   integer, allocatable :: sweepII_i(:), sweepII_j(:)


   !Adding sweep data
   INTEGER, PARAMETER :: total_max = 1000
   INTEGER, DIMENSION(1:total_max) :: sweep_track = 0
   INTEGER :: sweep_count = 0, total_sweep = 0, max_sweep = 0
   real(8):: avg_sweep

   INTEGER, DIMENSION(1:total_max) :: neg_track = 0
   INTEGER, DIMENSION(1:total_max) :: pos_track = 0
   INTEGER :: t_neg_count, t_pos_count, total_neg, total_pos, max_neg, max_pos
   real(8):: avg_pos, avg_neg



   if (mpirank == 0) then
      allocate (sweepI_i(total_x*total_y), sweepI_j(total_x*total_y))
      allocate (sweepII_i(total_x*total_y), sweepII_j(total_x*total_y))
      call sweeping_order
   end if

   !if (mpirank == 0) then
   !   allocate (u_global(4, -2:total_x + 3, -2:total_y + 3))
   !end if

   DO j = 1, total_y
      DO i = 1, total_x
         weights(i, j) = 1./(nSteps_x*nSteps_y)
      END DO
   END DO

   call get_command_argument(1, output_dir)

   call MPI_Init(mpierror)
   call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpierror)
   call MPI_Comm_size(MPI_COMM_WORLD, mpisize, mpierror)

   IF (mpisize .ne. cut_x*cut_y) then
      write (*, *) 'Need ', cut_x*cut_y, ' nodes!!!'
      stop
   END IF

   call initdata

   time = 0d0
   kcount = 0
   time_start = MPI_Wtime()
   DO WHILE (time <= time_terminal - 1d-10)
      CALL rk
      time = time + dt
      kcount = kcount + 1
      if (mpirank == 0 .and. kcount/50*50 .eq. kcount) write (6, *) kcount, time/time_terminal*100d0, '%'
   END DO
   time_end = MPI_Wtime()
   if (mpirank == 0) print *, 'CPU Time', time_end - time_start, kcount, time
   CALL output
   call MPI_Finalize(mpierror)

contains
   !***********************************************
   FUNCTION fflux(a)
      IMPLICIT NONE
      REAL, INTENT(IN):: a(5)
      REAL:: vx, vy, pres, rho, mon_x, mon_y, e, y
      REAL, DIMENSION(1:5)::fflux

      rho = a(1)
      mon_x = a(2)
      mon_y = a(3)
      e = a(4)
      y = a(5)

      vx = mon_x/rho
      vy = mon_y/rho

      if (rho < eps2) vx = 0.
      if (rho < eps2) vy = 0.

      pres = eos(rho, vx, vy, e, y)

      fflux(1) = mon_x
      fflux(2) = rho*vx**2 + pres
      fflux(3) = rho*vx*vy
      fflux(4) = (e + pres)*vx
      fflux(5) = y*vx

      RETURN
   END FUNCTION fflux

   FUNCTION gflux(a)
      IMPLICIT NONE
      REAL, INTENT(IN):: a(5)
      REAL:: vx, vy, pres, rho, mon_x, mon_y, e, y
      REAL, DIMENSION(1:5)::gflux

      rho = a(1)
      mon_x = a(2)
      mon_y = a(3)
      e = a(4)
      y = a(5)

      vx = mon_x/rho
      vy = mon_y/rho

      if (rho < eps2) vx = 0.
      if (rho < eps2) vy = 0.

      pres = eos(rho, vx, vy, e, y)

      gflux(1) = mon_y
      gflux(2) = rho*vx*vy
      gflux(3) = rho*vy**2 + pres
      gflux(4) = (e + pres)*vy
      gflux(5) = y*vy

      RETURN
   END FUNCTION gflux
   !****************************
   function eos(rho, vx, vy, e, y)
      implicit none
      REAL, INTENT(IN):: rho, vx, vy, e, y
      real:: eos

      eos = gamma1*(e - 0.5*rho*vx**2 - 0.5*rho*vy**2 - q_const*y)
   end function eos
   !**************************
   SUBROUTINE wenofd
      IMPLICIT NONE
      REAL:: evr(5, 5), evl(5, 5), cu_p(5, -2:2), cu_n(5, -2:2), ma(5, 5)
      REAL, DIMENSION(0:max(nx, ny) + 2):: w, vx, vy, h, ratio
      REAL:: den, xmt, ymt, eng, t0, t1, t2, t3, pre, vex, vey, cvel, vxm, vym, hm, qm, cm, rcm, b1, b2, ym, yrho
      REAL:: ul(5), ur(5), d_ave, m_ave, e_ave, p_ave, theta, point(5, 2), d_l, m_l, e_l, p_l, n_l, n_ave, y_l, y_ave
      real:: alpha_local

      call findalpha

      DO j = 0, ny + 1
         if (position_x == 0) then
            if (position_y*ny + j <= total*step) then
               do i = 0, 2
                  u(:, total - i, j) = u(:, total + 1 + i, j)
                  u(2, total - i, j) = -u(2, total + 1 + i, j)
               end do
            else
            end if
         else
         end if

         do i = 0, nx + 2
            den = u(1, i, j)
            xmt = u(2, i, j)
            ymt = u(3, i, j)
            eng = u(4, i, j)
            yrho = u(5, i, j)

            t0 = 1./den
            vex = xmt*t0
            vey = ymt*t0

            pre = gamma1*(eng - 0.5*(xmt*vex + ymt*vey) - q_const*yrho)
            w(i) = sqrt(abs(den))
            vx(i) = vex
            vy(i) = vey
            ratio(i) = yrho*t0
            h(i) = (pre + eng)*t0
         end do
         alpha_local = 0.
         do i = -2, nx + 3
            do ig = -mg, mg
               alpha_local = max(alpha_local, alpha(i, j + ig))
            end do
         end do

         do i = -2, nx + 3
            fp(:, i, j) = 0.5*(u(:, i, j) + fflux(u(:, i, j))/alpha_local)
            fn(:, i, j) = 0.5*(u(:, i, j) - fflux(u(:, i, j))/alpha_local)
         end do
         DO i = 0, nx + 1

            !********************************
            ! Roe Average Matrix for f'(u)
            !********************************

            t0 = w(i)/(w(i) + w(i + 1) + 1d-15)
            t1 = 1.-t0
            vxm = t0*vx(i) + t1*vx(i + 1)
            vym = t0*vy(i) + t1*vy(i + 1)
            ym = t0*ratio(i) + t1*ratio(i + 1)

            hm = t0*h(i) + t1*h(i + 1)
            qm = 0.5*(vxm*vxm + vym*vym)
            cm = sqrt(abs(gamma1*(hm - qm - q_const*ym))) + 1d-15
            t0 = vxm*cm
            evr(1, 1) = 1.0
            evr(1, 2) = 1.0
            evr(1, 3) = 0.0
            evr(1, 4) = 0.0
            evr(1, 5) = 1.0

            evr(2, 1) = vxm - cm
            evr(2, 2) = vxm
            evr(2, 3) = 0.
            evr(2, 4) = 0.
            evr(2, 5) = vxm + cm

            evr(3, 1) = vym
            evr(3, 2) = vym
            evr(3, 3) = -1.
            evr(3, 4) = 0.
            evr(3, 5) = vym

            evr(4, 1) = hm - t0
            evr(4, 2) = qm
            evr(4, 3) = -vym
            evr(4, 4) = q_const
            evr(4, 5) = hm + t0

            evr(5, 1) = ym
            evr(5, 2) = 0.
            evr(5, 3) = 0.
            evr(5, 4) = 1.
            evr(5, 5) = ym

            rcm = 1./cm
            b1 = gamma1*rcm*rcm
            b2 = qm*b1
            t0 = vxm*rcm
            t1 = b1*vxm
            t2 = 0.5*b1
            t3 = b1*vym
            evl(1, 1) = 0.5*(b2 + t0)
            evl(1, 2) = -0.5*(t1 + rcm)
            evl(1, 3) = -0.5*t3
            evl(1, 4) = t2
            evl(1, 5) = -t2*q_const

            evl(2, 1) = 1.-b2
            evl(2, 2) = t1
            evl(2, 3) = t3
            evl(2, 4) = -b1
            evl(2, 5) = b1*q_const

            evl(3, 1) = vym
            evl(3, 2) = 0.0
            evl(3, 3) = -1.0
            evl(3, 4) = 0.0
            evl(3, 5) = 0.0

            evl(4, 1) = -b2*ym
            evl(4, 2) = t1*ym
            evl(4, 3) = t3*ym
            evl(4, 4) = -b1*ym
            evl(4, 5) = 1.+q_const*ym*b1

            evl(5, 1) = 0.5*(b2 - t0)
            evl(5, 2) = -0.5*(t1 - rcm)
            evl(5, 3) = -0.5*t3
            evl(5, 4) = t2
            evl(5, 5) = -t2*q_const

            !******************************************
            ! projection into characteristic fields
            !******************************************
            DO ig = -2, 2
            DO k = 1, 5
               cu_p(k, ig) = evl(k, 1)*fp(1, i + ig, j) + evl(k, 2)*fp(2, i + ig, j) + evl(k, 3)*fp(3, i + ig, j)&
                         &+ evl(k, 4)*fp(4, i + ig, j) + evl(k, 5)*fp(5, i + ig, j)
               cu_n(k, ig) = evl(k, 1)*fn(1, i + ig, j) + evl(k, 2)*fn(2, i + ig, j) + evl(k, 3)*fn(3, i + ig, j)&
                         &+ evl(k, 4)*fn(4, i + ig, j) + evl(k, 5)*fn(5, i + ig, j)
            END DO
            END DO
            !**************************
            ! weno reconstruction
            !**************************
            DO kk = 1, 5
               call weno_right(cu_p(kk, :), ur(kk))
               call weno_left(cu_n(kk, :), ul(kk))
            END DO
            !******************************************
            ! projection back into conserved variables
            !******************************************
            DO k = 1, 5
               fp_flux(k, i, j) = evr(k, 1)*ur(1) + evr(k, 2)*ur(2) + evr(k, 3)*ur(3) + evr(k, 4)*ur(4) + evr(k, 5)*ur(5)
               fn_flux(k, i, j) = evr(k, 1)*ul(1) + evr(k, 2)*ul(2) + evr(k, 3)*ul(3) + evr(k, 4)*ul(4) + evr(k, 5)*ul(5)
            END DO

            if (use_limiter == 1) then

               !*****************************
               ! positivity limiter
               !*****************************
!          if (fp(1,i,j)<0. .or. fp(4,i,j)-fp(2,i,j)**2/fp(1,i,j)/2.-fp(3,i,j)**2/fp(1,i,j)/2.<0.) then
!          else
!          endif
!
!          if (fn(1,i,j)<0. .or. fn(4,i,j)-fn(2,i,j)**2/fn(1,i,j)/2.-fn(3,i,j)**2/fn(1,i,j)/2.<0.) then
!          else
!          endif

               theta = 1d0
               d_ave = fp(1, i, j)
               m_ave = fp(2, i, j)
               n_ave = fp(3, i, j)
               e_ave = fp(4, i, j)
               y_ave = fp(5, i, j)
               p_ave = gamma1*(e_ave - 0.5*m_ave**2/d_ave - 0.5*n_ave**2/d_ave - q_const*y_ave)
               theta = 1.
               point(:, 1) = fp_flux(:, i, j)
               point(:, 2) = (fp(:, i, j) - fp_flux(:, i, j)/12.)/(11./12.)

               do ig = 1, 2
                  d_l = point(1, ig)
                  if (d_l < 0.) theta = min(theta, abs(d_ave - 0.)/abs(d_ave - d_l))
               end do
               theta = max(theta, 0.)

               if (theta < 1.) then
                  fp_flux(1, i, j) = theta*fp_flux(1, i, j) + (1.-theta)*d_ave
                  den_count = den_count + 1
               else
               end if

               theta = 1.
               point(1, 1) = fp_flux(1, i, j)
               point(1, 2) = (fp(1, i, j) - fp_flux(1, i, j)/12.)/(11./12.)

               do ig = 1, 2
                  d_l = point(1, ig)
                  m_l = point(2, ig)
                  n_l = point(3, ig)
                  e_l = point(4, ig)
                  y_l = point(5, ig)

                  p_l = gamma1*(e_l - m_l**2/d_l/2.-n_l**2/d_l/2.-q_const*y_l)
                  if (p_l < 0.) theta = min(theta, abs(p_ave - 0.)/abs(p_ave - p_l))
               end do

               theta = max(theta, 0.)
               ! if (u(1,i,j)<1d-3) theta=0.
               if (theta < 1d0) then
                  pre_count = pre_count + 1
                  fp_flux(1, i, j) = theta*fp_flux(1, i, j) + (1.-theta)*d_ave
                  fp_flux(2, i, j) = theta*fp_flux(2, i, j) + (1.-theta)*m_ave
                  fp_flux(3, i, j) = theta*fp_flux(3, i, j) + (1.-theta)*n_ave
                  fp_flux(4, i, j) = theta*fp_flux(4, i, j) + (1.-theta)*e_ave
                  fp_flux(5, i, j) = theta*fp_flux(5, i, j) + (1.-theta)*y_ave

               else
               end if

               theta = 1d0
               d_ave = fn(1, i, j)
               m_ave = fn(2, i, j)
               n_ave = fn(3, i, j)
               e_ave = fn(4, i, j)
               y_ave = fn(5, i, j)
               p_ave = gamma1*(e_ave - 0.5*m_ave**2/d_ave - 0.5*n_ave**2/d_ave - q_const*y_ave)
               theta = 1.
               point(:, 1) = fn_flux(:, i, j)
               point(:, 2) = (fn(:, i, j) - fn_flux(:, i, j)/12.)/(11./12.)

               do ig = 1, 2
                  d_l = point(1, ig)
                  if (d_l < 0.) theta = min(theta, abs(d_ave - 0.)/abs(d_ave - d_l))
               end do
               theta = max(theta, 0.)

               if (theta < 1.) then
                  fn_flux(1, i, j) = theta*fn_flux(1, i, j) + (1.-theta)*d_ave
                  den_count = den_count + 1
               else
               end if

               theta = 1.
               point(1, 1) = fn_flux(1, i, j)
               point(1, 2) = (fn(1, i, j) - fn_flux(1, i, j)/12.)/(11./12.)

               do ig = 1, 2
                  d_l = point(1, ig)
                  m_l = point(2, ig)
                  n_l = point(3, ig)
                  e_l = point(4, ig)
                  y_l = point(5, ig)

                  p_l = gamma1*(e_l - m_l**2/d_l/2.-n_l**2/d_l/2.-q_const*y_l)
                  if (p_l < 0.) theta = min(theta, abs(p_ave - 0.)/abs(p_ave - p_l))
               end do

               theta = max(theta, 0.)
               ! if (u(1,i,j)<1d-3) theta=0.
               if (theta < 1d0) then
                  pre_count = pre_count + 1
                  fn_flux(1, i, j) = theta*fn_flux(1, i, j) + (1.-theta)*d_ave
                  fn_flux(2, i, j) = theta*fn_flux(2, i, j) + (1.-theta)*m_ave
                  fn_flux(3, i, j) = theta*fn_flux(3, i, j) + (1.-theta)*n_ave
                  fn_flux(4, i, j) = theta*fn_flux(4, i, j) + (1.-theta)*e_ave
                  fn_flux(5, i, j) = theta*fn_flux(5, i, j) + (1.-theta)*y_ave

               else
               end if
            end if

            fp_flux(:, i, j) = alpha_local*fp_flux(:, i, j)
            fn_flux(:, i, j) = alpha_local*fn_flux(:, i, j)
         END DO
      END DO

      if (position_x == 0) then
         do j = -2, ny + 3
            if (position_y*ny + j == total*step) then
               do i = 0, 2
               do k = 1, total
                  u(:, k, j - i) = u(:, k, j + 1 + i)
                  u(3, k, j - i) = -u(3, k, j + 1 + i)
               end do
               end do
            else
            end if

         end do
      else
      end if

      DO i = 0, nx + 1

         do j = 0, ny + 2
            den = u(1, i, j)
            xmt = u(2, i, j)
            ymt = u(3, i, j)
            eng = u(4, i, j)
            yrho = u(5, i, j)

            t0 = 1./den
            vex = xmt*t0
            vey = ymt*t0

            pre = gamma1*(eng - 0.5*(xmt*vex + ymt*vey) - q_const*yrho)
            w(j) = sqrt(abs(den))
            vx(j) = vex
            vy(j) = vey
            ratio(j) = yrho*t0
            h(j) = (pre + eng)*t0
         end do

         alpha_local = 0.
         do j = -2, ny + 3
            do ig = -mg, mg
               alpha_local = max(alpha_local, alpha(i + ig, j))
            end do
         end do
         do j = -2, ny + 3
            gp(:, i, j) = 0.5*(u(:, i, j) + gflux(u(:, i, j))/alpha_local)
            gn(:, i, j) = 0.5*(u(:, i, j) - gflux(u(:, i, j))/alpha_local)
         end do

         DO j = 0, ny + 1

            !********************************
            ! Roe Average Matrix for g'(u)
            !********************************

            t0 = w(j)/(w(j) + w(j + 1) + 1d-15)
            t1 = 1.-t0
            vxm = t0*vx(j) + t1*vx(j + 1)
            vym = t0*vy(j) + t1*vy(j + 1)
            hm = t0*h(j) + t1*h(j + 1)
            ym = t0*ratio(j) + t1*ratio(j + 1)

            qm = 0.5*(vxm*vxm + vym*vym)
            cm = sqrt(abs(gamma1*(hm - qm - q_const*ym))) + 1d-15
            t0 = vym*cm
            evr(1, 1) = 1.0
            evr(1, 2) = 1.0
            evr(1, 3) = 0.0
            evr(1, 4) = 0.0
            evr(1, 5) = 1.0

            evr(2, 1) = vxm
            evr(2, 2) = vxm
            evr(2, 3) = 1.0
            evr(2, 4) = 0.
            evr(2, 5) = vxm

            evr(3, 1) = vym - cm
            evr(3, 2) = vym
            evr(3, 3) = 0.
            evr(3, 4) = 0.
            evr(3, 5) = vym + cm

            evr(4, 1) = hm - t0
            evr(4, 2) = qm
            evr(4, 3) = vxm
            evr(4, 4) = q_const
            evr(4, 5) = hm + t0

            evr(5, 1) = ym
            evr(5, 2) = 0.
            evr(5, 3) = 0.
            evr(5, 4) = 1.
            evr(5, 5) = ym

            rcm = 1./cm
            b1 = gamma1*rcm*rcm
            b2 = qm*b1
            t0 = vym*rcm
            t1 = b1*vym
            t2 = 0.5*b1
            t3 = b1*vxm
            evl(1, 1) = 0.5*(b2 + t0)
            evl(1, 2) = -0.5*t3
            evl(1, 3) = -0.5*(t1 + rcm)
            evl(1, 4) = t2
            evl(1, 5) = -q_const*t2

            evl(2, 1) = 1.-b2
            evl(2, 2) = t3
            evl(2, 3) = t1
            evl(2, 4) = -b1
            evl(2, 5) = b1*q_const

            evl(3, 1) = -vxm
            evl(3, 2) = 1.0
            evl(3, 3) = 0.0
            evl(3, 4) = 0.0
            evl(3, 5) = 0.0

            evl(4, 1) = -ym*b2
            evl(4, 2) = ym*t3
            evl(4, 3) = ym*t1
            evl(4, 4) = -ym*b1
            evl(4, 5) = 1.+ym*b1*q_const

            evl(5, 1) = 0.5*(b2 - t0)
            evl(5, 2) = -0.5*t3
            evl(5, 3) = -0.5*(t1 - rcm)
            evl(5, 4) = t2
            evl(5, 5) = -q_const*t2

            !******************************************
            ! projection into characteristic fields
            !******************************************
            DO ig = -2, 2
            DO k = 1, 5
               cu_p(k, ig) = evl(k, 1)*gp(1, i, j + ig) + evl(k, 2)*gp(2, i, j + ig) + evl(k, 3)*gp(3, i, j + ig)&
                         &+ evl(k, 4)*gp(4, i, j + ig) + evl(k, 5)*gp(5, i, j + ig)
               cu_n(k, ig) = evl(k, 1)*gn(1, i, j + ig) + evl(k, 2)*gn(2, i, j + ig) + evl(k, 3)*gn(3, i, j + ig)&
                         &+ evl(k, 4)*gn(4, i, j + ig) + evl(k, 5)*gn(5, i, j + ig)
            END DO
            END DO
            !**************************
            ! weno reconstruction
            !**************************
            DO kk = 1, 5
               call weno_right(cu_p(kk, :), ur(kk))
               call weno_left(cu_n(kk, :), ul(kk))
            END DO
            !******************************************
            ! projection back into conserved variables
            !******************************************
            DO k = 1, 5
               gp_flux(k, i, j) = evr(k, 1)*ur(1) + evr(k, 2)*ur(2) + evr(k, 3)*ur(3) + evr(k, 4)*ur(4) + evr(k, 5)*ur(5)
               gn_flux(k, i, j) = evr(k, 1)*ul(1) + evr(k, 2)*ul(2) + evr(k, 3)*ul(3) + evr(k, 4)*ul(4) + evr(k, 5)*ul(5)
            END DO

            if (use_limiter == 1) then
               !*****************************
               ! positivity limiter
               !*****************************
!          if (gp(1,i,j)<0. .or. gp(4,i,j)-gp(2,i,j)**2/gp(1,i,j)/2.-gp(3,i,j)**2/gp(1,i,j)/2.<0.) then
!          else
!          endif
!
!          if (gn(1,i,j)<0. .or. gn(4,i,j)-gn(2,i,j)**2/gn(1,i,j)/2.-gn(3,i,j)**2/gn(1,i,j)/2.<0.) then
!          else
!          endif

               theta = 1d0
               d_ave = gp(1, i, j)
               m_ave = gp(2, i, j)
               n_ave = gp(3, i, j)
               e_ave = gp(4, i, j)
               y_ave = gp(5, i, j)
               p_ave = gamma1*(e_ave - 0.5*m_ave**2/d_ave - 0.5*n_ave**2/d_ave - q_const*y_ave)
               theta = 1.
               point(:, 1) = gp_flux(:, i, j)
               point(:, 2) = (gp(:, i, j) - gp_flux(:, i, j)/12.)/(11./12.)

               do ig = 1, 2
                  d_l = point(1, ig)
                  if (d_l < 0.) theta = min(theta, abs(d_ave - 0.)/abs(d_ave - d_l))
               end do
               theta = max(theta, 0.)

               if (theta < 1.) then
                  gp_flux(1, i, j) = theta*gp_flux(1, i, j) + (1.-theta)*d_ave
                  den_count = den_count + 1
               else
               end if

               theta = 1.
               point(1, 1) = gp_flux(1, i, j)
               point(1, 2) = (gp(1, i, j) - gp_flux(1, i, j)/12.)/(11./12.)

               do ig = 1, 2
                  d_l = point(1, ig)
                  m_l = point(2, ig)
                  n_l = point(3, ig)
                  e_l = point(4, ig)
                  y_l = point(5, ig)

                  p_l = gamma1*(e_l - m_l**2/d_l/2.-n_l**2/d_l/2.-q_const*y_l)
                  if (p_l < 0.) theta = min(theta, abs(p_ave - 0.)/abs(p_ave - p_l))
               end do

               theta = max(theta, 0.)
               ! if (u(1,i,j)<1d-3) theta=0.
               if (theta < 1d0) then
                  pre_count = pre_count + 1
                  gp_flux(1, i, j) = theta*gp_flux(1, i, j) + (1.-theta)*d_ave
                  gp_flux(2, i, j) = theta*gp_flux(2, i, j) + (1.-theta)*m_ave
                  gp_flux(3, i, j) = theta*gp_flux(3, i, j) + (1.-theta)*n_ave
                  gp_flux(4, i, j) = theta*gp_flux(4, i, j) + (1.-theta)*e_ave
                  gp_flux(5, i, j) = theta*gp_flux(5, i, j) + (1.-theta)*y_ave

               else
               end if

               theta = 1d0
               d_ave = gn(1, i, j)
               m_ave = gn(2, i, j)
               n_ave = gn(3, i, j)
               e_ave = gn(4, i, j)
               y_ave = gn(5, i, j)
               p_ave = gamma1*(e_ave - 0.5*m_ave**2/d_ave - 0.5*n_ave**2/d_ave - q_const*y_ave)
               theta = 1.
               point(:, 1) = gn_flux(:, i, j)
               point(:, 2) = (gn(:, i, j) - gn_flux(:, i, j)/12.)/(11./12.)

               do ig = 1, 2
                  d_l = point(1, ig)
                  if (d_l < 0.) theta = min(theta, abs(d_ave - 0.)/abs(d_ave - d_l))
               end do
               theta = max(theta, 0.)

               if (theta < 1.) then
                  gn_flux(1, i, j) = theta*gn_flux(1, i, j) + (1.-theta)*d_ave
                  den_count = den_count + 1
               else
               end if

               theta = 1.
               point(1, 1) = gn_flux(1, i, j)
               point(1, 2) = (gn(1, i, j) - gn_flux(1, i, j)/12.)/(11./12.)

               do ig = 1, 2
                  d_l = point(1, ig)
                  m_l = point(2, ig)
                  n_l = point(3, ig)
                  e_l = point(4, ig)
                  y_l = point(5, ig)

                  p_l = gamma1*(e_l - m_l**2/d_l/2.-n_l**2/d_l/2.-q_const*y_l)
                  if (p_l < 0.) theta = min(theta, abs(p_ave - 0.)/abs(p_ave - p_l))
               end do

               theta = max(theta, 0.)
               ! if (u(1,i,j)<1d-3) theta=0.
               if (theta < 1d0) then
                  pre_count = pre_count + 1
                  gn_flux(1, i, j) = theta*gn_flux(1, i, j) + (1.-theta)*d_ave
                  gn_flux(2, i, j) = theta*gn_flux(2, i, j) + (1.-theta)*m_ave
                  gn_flux(3, i, j) = theta*gn_flux(3, i, j) + (1.-theta)*n_ave
                  gn_flux(4, i, j) = theta*gn_flux(4, i, j) + (1.-theta)*e_ave
                  gn_flux(5, i, j) = theta*gn_flux(5, i, j) + (1.-theta)*y_ave

               else
               end if

            end if

            gp_flux(:, i, j) = alpha_local*gp_flux(:, i, j)
            gn_flux(:, i, j) = alpha_local*gn_flux(:, i, j)

         END DO
      END DO

   END SUBROUTINE wenofd
   !******************************
   SUBROUTINE scheme
      IMPLICIT NONE

      DO i = 1, nx
      DO j = 1, ny
         RHS(:, i, j) = -((fp_flux(:, i, j) - fn_flux(:, i + 1, j)) - (fp_flux(:, i - 1, j) - fn_flux(:, i, j)))/dx&
         &        - ((gp_flux(:, i, j) - gn_flux(:, i, j + 1)) - (gp_flux(:, i, j - 1) - gn_flux(:, i, j)))/dy
         RHS(5, i, j) = RHS(5, i, j) + chemsource(u(:, i, j))

      END DO
      END DO

   END SUBROUTINE scheme
   !**********************************
   function chemsource(a)
      implicit none
      REAL, INTENT(IN):: a(5)
      real:: chemsource, vx, vy, rho, mon_x, mon_y, e, y, pres

      rho = a(1)
      mon_x = a(2)
      mon_y = a(3)
      e = a(4)
      y = a(5)

      vx = mon_x/rho
      vy = mon_y/rho

      if (rho < 1d-8) vx = 0.
      if (rho < 1d-8) vy = 0.

      pres = eos(rho, vx, vy, e, y)
      chemsource = -k_const*y*exp(-t_const*rho/pres)
   end function chemsource
   !**********************

   SUBROUTINE findalpha
      IMPLICIT NONE
      REAL::den, pres, sound, alphatmp1, v_x, v_y, source_tmp

      source_alpha = 0.
      alphamax = 0.
      alpha = 0.
      do i = -2, nx + 3
      do j = -2, ny + 3
         den = u(1, i, j)
         if (den >= eps2) then
            v_x = u(2, i, j)/den
            v_y = u(3, i, j)/den
         else
            v_x = 0.
            v_y = 0.
         end if

         pres = gamma1*(u(4, i, j) - 0.5*den*v_x**2 - 0.5*den*v_y**2 - u(5, i, j)*q_const)
         sound = sqrt(gamma*pres/den)

         alpha(i, j) = max(abs(u(2, i, j)/u(1, i, j)) + sound, abs(u(3, i, j)/u(1, i, j)) + sound)
         alphamax = max(alphamax, alpha(i, j))
         source_alpha = max(source_alpha, k_const*exp(-t_const*den/pres))
      end do
      end do

      call MPI_AllReduce(alphamax, alphatmp1, 1, MPI_REAL8, MPI_MAX, &
      & MPI_COMM_WORLD, mpierror)
      alphamax = alphatmp1
      call MPI_AllReduce(source_alpha, source_tmp, 1, MPI_REAL8, MPI_MAX, &
      & MPI_COMM_WORLD, mpierror)
      source_alpha = source_tmp

   END SUBROUTINE findalpha
   !**************************
   SUBROUTINE rk
      IMPLICIT NONE
      REAL(8), PARAMETER:: cfl = 0.88
      REAL(8):: alphatemp, dt_cfl_x, dt_cfl_y, dt_cfl
      rktag = 1
      CALL wenofd

      dt_cfl_x = dx*cfl/alphamax
      dt_cfl_y = dy*cfl/alphamax
      dt_cfl = min(dt_cfl_x, dt_cfl_y)

      dt = min(dt_cfl, time_terminal - time)

      alphatemp = alphamax
      postag = 0

      rktag = 1
      u0 = u
      CALL scheme
      u = u + dt*RHS

      call checkpositivity

      call gather_to_root

      if (use_limiter == 2) then
         if (postag == 1) then
            if (mpirank == 0) then
               call sweep_lim
            end if
         end if
      end if

      call scatter_from_root

      CALL bc

      rktag = 2
      CALL wenofd

      CALL scheme
      u = 0.75d0*u0 + 0.25d0*(u + dt*RHS)
      call checkpositivity

      call gather_to_root
      if (use_limiter == 2) then
         if (postag == 1) then
            if (mpirank == 0) then
               call sweep_lim
            end if
         end if
      end if

      call scatter_from_root

      CALL bc

      rktag = 3
      CALL wenofd

      CALL scheme
      u = (u0 + 2d0*(u + dt*RHS))/3d0

      call checkpositivity

      call gather_to_root
      if (use_limiter == 2) then
         if (postag == 1) then
            if (mpirank == 0) then
               call sweep_lim
            end if
         end if
      end if

      call scatter_from_root

      CALL bc

   END SUBROUTINE rk
   !*************************
   SUBROUTINE checkpositivity
      IMPLICIT NONE
      integer:: postmp
      postag = 0
      do i = 1, nx
      do j = 1, ny
         if (u(1, i, j) <= eps) then
            postag = 1
         else
         end if

         if (u(4, i, j) - u(3, i, j)**2/u(1, i, j)/2.-u(2, i, j)**2/u(1, i, j)/2.-q_const*u(5, i, j) <= eps) then
            postag = 1
         else
         end if
      end do
      end do


      call MPI_AllReduce(postag, postmp, 1, MPI_INTEGER, MPI_MAX, &
      & MPI_COMM_WORLD, mpierror)

      postag = postmp
   END SUBROUTINE checkpositivity
   !*************************
   SUBROUTINE bc
      IMPLICIT NONE
      integer, parameter:: total_elems = 2*(15*(ny + 6)) + 2*(15*(nx + 6))
      integer, parameter:: commsize = 15*(ny + 2 + nx + 2 + 6)*100
      real, dimension(commsize) :: buffer
      real, dimension(5, 3, -2:ny + 3):: send_buf_y, recv_buf_y
      real, dimension(5, -2:nx + 3, 3):: send_buf_x, recv_buf_x
      INTEGER :: buf_size
      buf_size = commsize*8

      if (cut_x == 1) then

      else
         call exchange_x
      end if
      if (cut_y == 1) then

      else
         call exchange_y
      end if

      !*************************************************
      ! do it twice to update the info in the 4 corners
      !*************************************************
      if (cut_x == 1) then

      else
         call exchange_x
      end if
      if (cut_y == 1) then

      else
         call exchange_y
      end if

      if (cut_x == 1 .and. cut_y == 1) then
         !call MPI_BUFFER_ATTACH(buffer, commsize, mpierror)
         call MPI_BUFFER_ATTACH(buffer, buf_size, mpierror)

         mpitag = 1
         send_buf_y(:, 1, :) = u(:, 1, :)
         send_buf_y(:, 2, :) = u(:, 2, :)
         send_buf_y(:, 3, :) = u(:, 3, :)
         call MPI_Bsend(send_buf_y, 15*(ny + 6), MPI_REAL8,&
        & destination(1), mpitag, MPI_COMM_WORLD, mpierror)
         mpitag = 2
         send_buf_y(:, 1, :) = u(:, nx, :)
         send_buf_y(:, 2, :) = u(:, nx - 1, :)
         send_buf_y(:, 3, :) = u(:, nx - 2, :)
         call MPI_Bsend(send_buf_y, 15*(ny + 6), MPI_REAL8,&
        & destination(2), mpitag, MPI_COMM_WORLD, mpierror)
         mpitag = 3
         send_buf_x(:, :, 1) = u(:, :, ny)
         send_buf_x(:, :, 2) = u(:, :, ny - 1)
         send_buf_x(:, :, 3) = u(:, :, ny - 2)
         call MPI_Bsend(send_buf_x, 15*(nx + 6), MPI_REAL8,&
        & destination(3), mpitag, MPI_COMM_WORLD, mpierror)
         mpitag = 4
         send_buf_x(:, :, 1) = u(:, :, 1)
         send_buf_x(:, :, 2) = u(:, :, 2)
         send_buf_x(:, :, 3) = u(:, :, 3)
         call MPI_Bsend(send_buf_x, 15*(nx + 6), MPI_REAL8,&
        & destination(4), mpitag, MPI_COMM_WORLD, mpierror)

         mpitag = 1
         call MPI_RECV(recv_buf_y, 15*(ny + 6), MPI_REAL8, &
        & source(2), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, nx + 1, :) = recv_buf_y(:, 1, :)
         u(:, nx + 2, :) = recv_buf_y(:, 2, :)
         u(:, nx + 3, :) = recv_buf_y(:, 3, :)
         mpitag = 2
         call MPI_RECV(recv_buf_y, 15*(ny + 6), MPI_REAL8, &
        & source(1), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, 0, :) = recv_buf_y(:, 1, :)
         u(:, -1, :) = recv_buf_y(:, 2, :)
         u(:, -2, :) = recv_buf_y(:, 3, :)

         mpitag = 3
         call MPI_RECV(recv_buf_x, 15*(nx + 6), MPI_REAL8, &
        & source(4), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, :, 0) = recv_buf_x(:, :, 1)
         u(:, :, -1) = recv_buf_x(:, :, 2)
         u(:, :, -2) = recv_buf_x(:, :, 3)
         mpitag = 4
         call MPI_RECV(recv_buf_x, 15*(nx + 6), MPI_REAL8, &
        & source(3), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, :, ny + 1) = recv_buf_x(:, :, 1)
         u(:, :, ny + 2) = recv_buf_x(:, :, 2)
         u(:, :, ny + 3) = recv_buf_x(:, :, 3)

         !call MPI_BUFFER_DETACH(buffer, commsize, mpierror)
         call MPI_BUFFER_DETACH(buffer, buf_size, mpierror)

      else
      end if

      !************************************************************
      ! outflow for other boundaries
      !************************************************************

      if (position_y == 0) then
         do i = -2, 0
            u(1, :, i) = 1.
            u(2, :, i) = 0.
            u(3, :, i) = 0.
            u(4, :, i) = 55.
            u(5, :, i) = 1.
         end do
      else
      end if
      if (position_x == cut_x - 1) then
         do i = 1, 3
            u(1, nx + i, :) = 1.
            u(2, nx + i, :) = 0.
            u(3, nx + i, :) = 0.
            u(4, nx + i, :) = 55.
            u(5, nx + i, :) = 1.
         end do

      else
      end if

      if (position_y == cut_y - 1) then
         do i = 1, 3
            u(:, :, ny + i) = u(:, :, ny)
            u(3, :, ny + i) = -u(3, :, ny)
         end do
      else
      end if
      if (position_x == 0) then
         do i = -2, 0

            u(1, i, :) = 11.
            u(2, i, :) = 11.*6.18
            u(3, i, :) = 0.
            u(4, i, :) = 970.
            u(5, i, :) = 11.

         end do

      else
      end if

      do i = -2, nx + 3
         do j = -2, ny + 3
            if (i + position_x*nx < total .and. j + position_y*ny < step*total) then
               u(1, i, j) = 1.
               u(2, i, j) = 0.
               u(3, i, j) = 0.
               u(4, i, j) = 55.
               u(5, i, j) = 1.
            else
            end if
         end do
      end do

   END SUBROUTINE bc

   SUBROUTINE exchange_x
      IMPLICIT NONE
      real, dimension(5, 3, -2:ny + 3):: send_buf_y, recv_buf_y

      IF (mod(position_x, 2) == 0) THEN
         mpitag = 1
         send_buf_y(:, 1, :) = u(:, 1, :)
         send_buf_y(:, 2, :) = u(:, 2, :)
         send_buf_y(:, 3, :) = u(:, 3, :)
         call MPI_send(send_buf_y, 15*(ny + 6), MPI_REAL8,&
        & destination(1), mpitag, MPI_COMM_WORLD, mpierror)

         mpitag = 1
         call MPI_RECV(recv_buf_y, 15*(ny + 6), MPI_REAL8, &
        & source(2), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, nx + 1, :) = recv_buf_y(:, 1, :)
         u(:, nx + 2, :) = recv_buf_y(:, 2, :)
         u(:, nx + 3, :) = recv_buf_y(:, 3, :)
      ELSE
         mpitag = 1
         call MPI_RECV(recv_buf_y, 15*(ny + 6), MPI_REAL8, &
        & source(2), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, nx + 1, :) = recv_buf_y(:, 1, :)
         u(:, nx + 2, :) = recv_buf_y(:, 2, :)
         u(:, nx + 3, :) = recv_buf_y(:, 3, :)
         mpitag = 1
         send_buf_y(:, 1, :) = u(:, 1, :)
         send_buf_y(:, 2, :) = u(:, 2, :)
         send_buf_y(:, 3, :) = u(:, 3, :)
         call MPI_send(send_buf_y, 15*(ny + 6), MPI_REAL8,&
        & destination(1), mpitag, MPI_COMM_WORLD, mpierror)
      END IF

      IF (mod(position_x, 2) == 0) THEN
         mpitag = 2
         send_buf_y(:, 1, :) = u(:, nx, :)
         send_buf_y(:, 2, :) = u(:, nx - 1, :)
         send_buf_y(:, 3, :) = u(:, nx - 2, :)
         call MPI_send(send_buf_y, 15*(ny + 6), MPI_REAL8,&
        & destination(2), mpitag, MPI_COMM_WORLD, mpierror)
         mpitag = 2
         call MPI_RECV(recv_buf_y, 15*(ny + 6), MPI_REAL8, &
        & source(1), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, 0, :) = recv_buf_y(:, 1, :)
         u(:, -1, :) = recv_buf_y(:, 2, :)
         u(:, -2, :) = recv_buf_y(:, 3, :)
      ELSE
         mpitag = 2
         call MPI_RECV(recv_buf_y, 15*(ny + 6), MPI_REAL8, &
        & source(1), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, 0, :) = recv_buf_y(:, 1, :)
         u(:, -1, :) = recv_buf_y(:, 2, :)
         u(:, -2, :) = recv_buf_y(:, 3, :)
         mpitag = 2
         send_buf_y(:, 1, :) = u(:, nx, :)
         send_buf_y(:, 2, :) = u(:, nx - 1, :)
         send_buf_y(:, 3, :) = u(:, nx - 2, :)
         call MPI_send(send_buf_y, 15*(ny + 6), MPI_REAL8,&
        & destination(2), mpitag, MPI_COMM_WORLD, mpierror)
      END IF
   END SUBROUTINE exchange_x

   SUBROUTINE exchange_y
      IMPLICIT NONE
      real, dimension(5, -2:nx + 3, 3):: send_buf_x, recv_buf_x

      IF (mod(position_y, 2) == 0) THEN
         mpitag = 4
         send_buf_x(:, :, 1) = u(:, :, 1)
         send_buf_x(:, :, 2) = u(:, :, 2)
         send_buf_x(:, :, 3) = u(:, :, 3)
         call MPI_send(send_buf_x, 15*(nx + 6), MPI_REAL8,&
        & destination(4), mpitag, MPI_COMM_WORLD, mpierror)
         mpitag = 4
         call MPI_RECV(recv_buf_x, 15*(nx + 6), MPI_REAL8, &
        & source(3), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, :, ny + 1) = recv_buf_x(:, :, 1)
         u(:, :, ny + 2) = recv_buf_x(:, :, 2)
         u(:, :, ny + 3) = recv_buf_x(:, :, 3)
      ELSE
         mpitag = 4
         call MPI_RECV(recv_buf_x, 15*(nx + 6), MPI_REAL8, &
        & source(3), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, :, ny + 1) = recv_buf_x(:, :, 1)
         u(:, :, ny + 2) = recv_buf_x(:, :, 2)
         u(:, :, ny + 3) = recv_buf_x(:, :, 3)
         mpitag = 4
         send_buf_x(:, :, 1) = u(:, :, 1)
         send_buf_x(:, :, 2) = u(:, :, 2)
         send_buf_x(:, :, 3) = u(:, :, 3)
         call MPI_send(send_buf_x, 15*(nx + 6), MPI_REAL8,&
        & destination(4), mpitag, MPI_COMM_WORLD, mpierror)
      END IF

      IF (mod(position_y, 2) == 0) THEN
         mpitag = 3
         send_buf_x(:, :, 1) = u(:, :, ny)
         send_buf_x(:, :, 2) = u(:, :, ny - 1)
         send_buf_x(:, :, 3) = u(:, :, ny - 2)
         call MPI_send(send_buf_x, 15*(nx + 6), MPI_REAL8,&
        & destination(3), mpitag, MPI_COMM_WORLD, mpierror)
         mpitag = 3
         call MPI_RECV(recv_buf_x, 15*(nx + 6), MPI_REAL8, &
        & source(4), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, :, 0) = recv_buf_x(:, :, 1)
         u(:, :, -1) = recv_buf_x(:, :, 2)
         u(:, :, -2) = recv_buf_x(:, :, 3)
      ELSE
         mpitag = 3
         call MPI_RECV(recv_buf_x, 15*(nx + 6), MPI_REAL8, &
        & source(4), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, :, 0) = recv_buf_x(:, :, 1)
         u(:, :, -1) = recv_buf_x(:, :, 2)
         u(:, :, -2) = recv_buf_x(:, :, 3)
         mpitag = 3
         send_buf_x(:, :, 1) = u(:, :, ny)
         send_buf_x(:, :, 2) = u(:, :, ny - 1)
         send_buf_x(:, :, 3) = u(:, :, ny - 2)
         call MPI_send(send_buf_x, 15*(nx + 6), MPI_REAL8,&
        & destination(3), mpitag, MPI_COMM_WORLD, mpierror)
      END IF
   END SUBROUTINE exchange_y
   !**************************
   subroutine initdata
      implicit none
      real:: xleft, xright, yleft, yright, ue(4)
      real rho_L, v_L, p_L, rho_R, v_R, p_R
      position_x = mod(mpirank, cut_x)
      position_y = mpirank/cut_x

      destination(1) = mpirank - 1
      destination(2) = mpirank + 1
      destination(3) = mpirank + cut_x
      destination(4) = mpirank - cut_x
!       destination(5)=mpirank-1+cut_x
!       destination(6)=mpirank+1+cut_x
!       destination(7)=mpirank-1-cut_x
!       destination(8)=mpirank+1-cut_x

      if (position_x == 0) destination(1) = mpirank - 1 + cut_x
      if (position_x == cut_x - 1) destination(2) = mpirank + 1 - cut_x
      if (position_y == 0) destination(4) = (cut_y - 1)*cut_x + mpirank
      if (position_y == cut_y - 1) destination(3) = mpirank - (cut_y - 1)*cut_x

!       if (position_x==0) destination(5)=mod(mpirank+2*cut_x-1, mpisize)
!       if (position_x==cut_x-1) destination(2)=mpirank+1-cut_x
!       if (position_y==0) destination(7)=mpirank-1
!       if (position_y==cut_y-1) destination(3)=mpirank-(cut_y-1)*cut_x

      source = destination

      xleft = 5./cut_x*position_x
      xright = 5./cut_x*(position_x + 1)
      xlen = xright - xleft
      dx = xlen/nx
      yleft = 5./cut_y*position_y
      yright = 5./cut_y*(position_y + 1)
      ylen = yright - yleft
      dy = ylen/ny

      do i = -2, nx + 3
         x(i) = xleft + i*dx - dx/2.
      end do
      do i = -2, ny + 3
         y(i) = yleft + i*dy - dy/2.
      end do

      u = 0.
      do i = 0, nx + 1
      do j = 0, ny + 1
         if (x(i) < 0.5) then
            u(1, i, j) = 11.
            u(2, i, j) = 11.*6.18
            u(3, i, j) = 0.
            u(4, i, j) = 970.
            u(5, i, j) = 11.
         else
            u(1, i, j) = 1.
            u(2, i, j) = 0.
            u(3, i, j) = 0.
            u(4, i, j) = 55.
            u(5, i, j) = 1.
         end if
      end do
      end do

      call bc
   end subroutine initdata
   !******************************

   !******************************
   SUBROUTINE output
      IMPLICIT NONE
      !real:: er1, eri, er1old, eriold, o1,oi, er1temp, eritemp, utemp(5,1:5*total,1:5*total), yrho

      real:: er1, eri, er1old, eriold, o1, oi, er1temp, eritemp, yrho

      REAL:: buffer(5*(cut_x*nx + 6)*(cut_y*ny + 6)), recv_buf(5, 1:nx, 1:ny), xtemp(1:total), ytemp(1:total)
      real:: lrho, lm, ln, le, pres
      integer:: den_tmp, pre_tmp
      integer :: buffer_size_bytes
      !Real(8), allocatable:: utemp(:,:,:)
      !ALLOCATE(utemp(5, 5*total, 5*total))

      !REAL(8):: recv_buf(5, 1:nx, 1:ny)
      !REAL(8), ALLOCATABLE :: buffer(:)
      !INTEGER(KIND=MPI_ADDRESS_KIND) :: buffer_addr
      !INTEGER :: buffer_size_bytes

      er1 = 0.
      eri = 0.

      buffer_size_bytes = (5*(cut_x*nx + 6)*(cut_y*ny + 6))*SIZEOF(1.0_8)

      !ALLOCATE(buffer(buffer_size_bytes))

      call MPI_AllReduce(den_count, den_tmp, 1, MPI_INTEGER, MPI_SUM, &
      & MPI_COMM_WORLD, mpierror)

      call MPI_AllReduce(pre_count, pre_tmp, 1, MPI_INTEGER, MPI_SUM, &
      & MPI_COMM_WORLD, mpierror)
      den_count = den_tmp
      pre_count = pre_tmp

      !call MPI_BUFFER_ATTACH(buffer, 10*(cut_x*nx+6)*(cut_y*ny+6), mpierror)
      !call MPI_BUFFER_ATTACH(buffer, buffer_size_bytes, mpierror)

      ! if (mpirank .ne.0) then
      !   mpitag=mpirank
      ! mpides=0
      !call MPI_Bsend(u(:,1:nx,1:ny),5*nx*ny,MPI_REAL8,mpides,&
      !        & mpitag,MPI_COMM_WORLD,mpierror)
      !  else
      !     do ig=1,nx
      !do jg=1,ny
      !        utemp(:,ig,jg)=u(:,ig,jg)
      !      enddo
      ! enddo
      ! DO i=1, mpisize-1
      !    mpisource=i
      !    mpitag=i
      !    call MPI_RECV(recv_buf, 5*nx*ny, &
      !    MPI_REAL8, mpisource, mpitag,MPI_COMM_WORLD, send_status,mpierror)
      !    DO ig=1,nx
      !    DO jg=1,ny
      !       utemp(:,ig+mod(i,cut_x)*nx,jg+i/cut_x*ny)=recv_buf(:,ig,jg)
      !    ENDDO
      !    ENDDO
      ! ENDDO
      !   endif
      !call MPI_BUFFER_DETACH(buffer, 10*(cut_x*nx+6)*(cut_y*ny+6), mpierror)
      !   call MPI_BUFFER_DETACH(buffer, buffer_size_bytes, mpierror)
      call gather_to_root

      if (mpirank == 0) then

         !open(73, file='den.dat')
         !open(74, file='pre.dat')
         !open(75, file='rate.dat')

         write (den_file_path, '(A,A)') trim(output_dir), "/den.dat"
         write (pre_file_path, '(A,A)') trim(output_dir), "/pre.dat"
         write (rate_file_path, '(A,A)') trim(output_dir), "/rate.dat"

         open (73, file=trim(den_file_path))
         open (74, file=trim(pre_file_path))
         open (75, file=trim(rate_file_path))

         write (73, *) 'zone i=', total*5, ', j=', total*5
         write (74, *) 'zone i=', total*5, ', j=', total*5
         write (75, *) 'zone i=', total*5, ', j=', total*5
         do j = 1, total*5
         do i = 1, total*5
            lrho = utemp(1, i, j)

            lm = utemp(2, i, j)
            ln = utemp(3, i, j)
            le = utemp(4, i, j)
            pres = gamma1*(le - 0.5*lm**2/lrho - 0.5*ln**2/lrho - q_const*utemp(5, i, j))
            yrho = utemp(5, i, j)/lrho
            if (i <= total .and. j < step*total) lrho = 100.
            if (i <= total .and. j < step*total) pres = 300.
            if (i <= total .and. j < step*total) yrho = 20.
            write (73, *) dx*(i - 0.5), dy*(j - 0.5), lrho
            write (74, *) dx*(i - 0.5), dy*(j - 0.5), pres
            write (75, *) dx*(i - 0.5), dy*(j - 0.5), yrho
         end do
         end do
         close (73)
         close (74)
         close (75)

         total_sweep = sum(sweep_track(1:sweep_count))
         avg_sweep = real(total_sweep,8)/real(sweep_count,8)
         max_sweep = maxval(sweep_track(1:sweep_count))


         print *, 'total sweep = ', total_sweep
         print *, 'avg sweep = ',   avg_sweep
         print *, 'max sweep = ',   max_sweep


         total_neg = sum(neg_track(1:sweep_count))
         avg_neg = real(total_neg,8)/real(sweep_count,8)
         max_neg = maxval(neg_track(1:sweep_count))

         print *, 'total negative pairs = ', total_neg
         print *, 'avg negative pairs = ',   avg_neg
         print *, 'max negative pairs = ',   max_neg

         total_pos = sum(pos_track(1:sweep_count))
         avg_pos = real(total_pos,8)/real(sweep_count,8)
         max_pos = maxval(pos_track(1:sweep_count))

         print *, 'total positive neighbors = ', total_pos
         print *, 'avg positive neighbors = ',   avg_pos
         print *, 'max positive neighbors = ',   max_pos
         
         print *, '# of limiter', den_count, pre_count
      else
      end if
   END SUBROUTINE output

      !Sweeping limiter routines added by Chloe Griffin

   !*****************************
   FUNCTION nod_lim(v_1, v_2, v_bar)
      IMPLICIT NONE
      REAL(8), DIMENSION(5), INTENT(IN):: v_1, v_2, v_bar
      REAL(8):: p_v_1, p_v_2, p_v_bar, t, t_1, delta, d_1_2, d_1_bar
      REAL(8), DIMENSION(5):: nod_lim

      P_v_1 = gamma1*(v_1(4) - 0.5*v_1(2)**2/v_1(1) - 0.5*v_1(3)**2/v_1(1) - q_const*v_1(5))
      P_v_2 = gamma1*(v_2(4) - 0.5*v_2(2)**2/v_2(1) - 0.5*v_2(3)**2/v_2(1) - q_const*v_2(5))

      !if (p_v_2 > 0 .AND. v_2(1)>0) then
      if (p_v_2 > eps) then
         t = ((P_v_1- eps)/(P_v_1 - P_v_2)); 
         t_pos_count = t_pos_count + 1 
      else
         P_v_bar = gamma1*(v_bar(4) - 0.5*v_bar(2)**2/v_bar(1) - 0.5*v_bar(3)**2/v_bar(1) - q_const*v_bar(5))
         ! t_1 = P_v_1/(P_v_1 - P_v_bar+1e-15);
         t_1 = (P_v_1-eps)/(P_v_1 - P_v_bar); 
         d_1_2 = sqrt((v_2(5) - v_1(5))**2 + (v_2(4) - v_1(4))**2 + (v_2(3) - v_1(3))**2 + (v_2(2) - v_1(2))**2 + (v_2(1) - v_1(1))**2); 
         d_1_bar = sqrt((v_bar(5) - v_1(5))**2 + (v_bar(4) - v_1(4))**2 + (v_bar(3) - v_1(3))**2 + (v_bar(2) - v_1(2))**2 + (v_bar(1) - v_1(1))**2); 
         t = min(t_1*d_1_bar/d_1_2, 0.25); 
         t_neg_count = t_neg_count + 1
      end if

      nod_lim(1) = (1 - t)*v_1(1) + t*v_2(1); 
      nod_lim(2) = (1 - t)*v_1(2) + t*v_2(2); 
      nod_lim(3) = (1 - t)*v_1(3) + t*v_2(3); 
      nod_lim(4) = (1 - t)*v_1(4) + t*v_2(4); 
      nod_lim(5) = (1 - t)*v_1(5) + t*v_2(5); 

      RETURN
   END FUNCTION nod_lim

   !*****************************
   !*****************************
   SUBROUTINE sweep_lim
      IMPLICIT NONE
      REAL(8), dimension(5)::  v_j, v_i, v_p, u_bar, v_j_t, v_i_t, v_p_t, d
      REAL(8), dimension(1:total_x, 1:total_y):: rho_t, m1_t, m2_t, E_t, P_check, ry_t
      REAL(8):: rho_avg, m1_avg, m2_avg, E_avg, ry_avg, pres_v_j, w_p, w_j, pres_v_i, w_i, rho_min, p_min
      INTEGER:: q

      REAL(8)::  pres_v_p, pres_v_i_t, pres_v_p_t, dist_tilde, dist_pi

!if (use_limiter == 2) then

      call scalar_sweep

      DO j = 1, total_y
         DO i = 1, total_x
            rho_t(i, j) = utemp(1, i, j); 
            m1_t(i, j) = utemp(2, i, j); 
            m2_t(i, j) = utemp(3, i, j); 
            E_t(i, j) = utemp(4, i, j); 
            ry_t(i, j) = utemp(5, i, j); 
         END DO
      END DO

      rho_avg = 0.0
      m1_avg = 0.0
      m2_avg = 0.0
      E_avg = 0.0
      ry_avg = 0.0

      do j = 1, total_y
         do i = 1, total_x
            rho_avg = rho_avg + weights(i, j)*rho_t(i, j)
            m1_avg = m1_avg + weights(i, j)*m1_t(i, j)
            m2_avg = m2_avg + weights(i, j)*m2_t(i, j)
            E_avg = E_avg + weights(i, j)*E_t(i, j)
            ry_avg = ry_avg + weights(i, j)*ry_t(i, j)
         end do
      end do

      u_bar(1) = rho_avg
      u_bar(2) = m1_avg
      u_bar(3) = m2_avg
      u_bar(4) = E_avg
      u_bar(5) = ry_avg

!gamma1*(u_bar(4) - 0.5*u_bar(2)**2/u_bar(1) - 0.5*u_bar(3)**2/u_bar(1) - q_const*u_bar(5))

!if (gamma1*(u_bar(4) - 0.5*(u_bar(2)*u_bar(2) + u_bar(3)*u_bar(3))/u_bar(1)) <= 0.) then
      if (gamma1*(E_avg - 0.5*m1_avg**2/rho_avg - 0.5*m2_avg**2/rho_avg - q_const*ry_avg) <= 0.) then
         print *, "ERROR: P(u_bar) <= 0."
         print *, "P(u_bar) = ", gamma1*(u_bar(4) - 0.5*u_bar(2)**2/u_bar(1) - 0.5*u_bar(3)**2/u_bar(1) - q_const*u_bar(5))
         stop "Program terminated due to incorrect physical conditions."
      end if

      sweep_count = sweep_count + 1
      t_pos_count = 0
      t_neg_count = 0

      DO q = 1, swp_num
         if (mod(q, 2) == 1) then  !Alternate sweep 1 and sweep 2
            !call scalar_sweep


            !Forward Sweep
            DO i = 1, (total_x*total_y - 1)

               v_i(1) = utemp(1, sweepI_i(i), sweepI_j(i)); 
               v_i(2) = utemp(2, sweepI_i(i), sweepI_j(i)); 
               v_i(3) = utemp(3, sweepI_i(i), sweepI_j(i)); 
               v_i(4) = utemp(4, sweepI_i(i), sweepI_j(i)); 
               v_i(5) = utemp(5, sweepI_i(i), sweepI_j(i)); 
               pres_v_i = gamma1*(v_i(4) - 0.5*v_i(2)*v_i(2)/v_i(1) - 0.5*v_i(3)**2/v_i(1) - q_const*v_i(5))
               w_i = weights(sweepI_i(i), sweepI_j(i)); 
               !if (pres_v_i <= 0. .OR. v_i(1) <= 0.) then
               if (pres_v_i <= eps .OR. v_i(1) <= eps) then

                  v_p(1) = utemp(1, sweepI_i(i + 1), sweepI_j(i + 1)); 
                  v_p(2) = utemp(2, sweepI_i(i + 1), sweepI_j(i + 1)); 
                  v_p(3) = utemp(3, sweepI_i(i + 1), sweepI_j(i + 1)); 
                  v_p(4) = utemp(4, sweepI_i(i + 1), sweepI_j(i + 1)); 
                  v_p(5) = utemp(5, sweepI_i(i + 1), sweepI_j(i + 1)); 
                  w_p = weights(sweepI_i(i + 1), sweepI_j(i + 1)); 
                  v_i_t = nod_lim(v_i, v_p, u_bar); 
                  d = v_i - v_i_t; 
                  v_p_t = v_p + (w_i/w_p)*(d); 
                  utemp(1, sweepI_i(i), sweepI_j(i)) = v_i_t(1); 
                  utemp(2, sweepI_i(i), sweepI_j(i)) = v_i_t(2); 
                  utemp(3, sweepI_i(i), sweepI_j(i)) = v_i_t(3); 
                  utemp(4, sweepI_i(i), sweepI_j(i)) = v_i_t(4); 
                  utemp(5, sweepI_i(i), sweepI_j(i)) = v_i_t(5); 
                  utemp(1, sweepI_i(i + 1), sweepI_j(i + 1)) = v_p_t(1); 
                  utemp(2, sweepI_i(i + 1), sweepI_j(i + 1)) = v_p_t(2); 
                  utemp(3, sweepI_i(i + 1), sweepI_j(i + 1)) = v_p_t(3); 
                  utemp(4, sweepI_i(i + 1), sweepI_j(i + 1)) = v_p_t(4); 
                  utemp(5, sweepI_i(i + 1), sweepI_j(i + 1)) = v_p_t(5); 
                  pres_v_p = gamma1*(v_p(4) - 0.5*v_p(2)*v_p(2)/v_p(1) - 0.5*v_p(3)**2/v_p(1) - q_const*v_p(5))
                  pres_v_i_t = gamma1*(v_i_t(4) - 0.5*v_i_t(2)**2/v_i_t(1) - 0.5*v_i_t(3)**2/v_i_t(1) - q_const*v_i_t(5))

                  pres_v_p_t = gamma1*(v_p_t(4) - 0.5*v_p_t(2)**2/v_p_t(1) - 0.5*v_p_t(3)**2/v_p_t(1) - q_const*v_p_t(5)); 
                  if (pres_v_i + pres_v_p > pres_v_i_t + pres_v_p_t + eps) then

                     ! Distance between v_p and v_i
                     dist_pi = sqrt(sum((v_p - v_i)**2))

                     ! Distance between v_p_t and v_i_t
                     dist_tilde = sqrt(sum((v_p_t - v_i_t)**2))

                     stop "Sweep check failed: total pressure increased after adjustment"
                  end if
                  ! if (i == 66 .AND. y_iter == 64) then
                  !end if
               end if
            end do

            !Backward Sweep
            DO i = total_x*total_y, 2, -1

               v_i(1) = utemp(1, sweepI_i(i), sweepI_j(i)); 
               v_i(2) = utemp(2, sweepI_i(i), sweepI_j(i)); 
               v_i(3) = utemp(3, sweepI_i(i), sweepI_j(i)); 
               v_i(4) = utemp(4, sweepI_i(i), sweepI_j(i)); 
               v_i(5) = utemp(5, sweepI_i(i), sweepI_j(i)); 
               pres_v_i = gamma1*(v_i(4) - 0.5*v_i(2)*v_i(2)/v_i(1) - 0.5*v_i(3)**2/v_i(1) - q_const*v_i(5))

               w_i = weights(sweepI_i(i), sweepI_j(i)); 
               !if (pres_v_i <= 0.  .OR. v_i(1) <= 0.) then
               if (pres_v_i <= eps .OR. v_i(1) <= eps) then

                  v_p(1) = utemp(1, sweepI_i(i - 1), sweepI_j(i - 1)); 
                  v_p(2) = utemp(2, sweepI_i(i - 1), sweepI_j(i - 1)); 
                  v_p(3) = utemp(3, sweepI_i(i - 1), sweepI_j(i - 1)); 
                  v_p(4) = utemp(4, sweepI_i(i - 1), sweepI_j(i - 1)); 
                  v_p(5) = utemp(5, sweepI_i(i - 1), sweepI_j(i - 1)); 
                  w_p = weights(sweepI_i(i - 1), sweepI_j(i - 1)); 
                  v_i_t = nod_lim(v_i, v_p, u_bar); 
                  d = v_i - v_i_t; 
                  v_p_t = v_p + (w_i/w_p)*(d); 
                  utemp(1, sweepI_i(i), sweepI_j(i)) = v_i_t(1); 
                  utemp(2, sweepI_i(i), sweepI_j(i)) = v_i_t(2); 
                  utemp(3, sweepI_i(i), sweepI_j(i)) = v_i_t(3); 
                  utemp(4, sweepI_i(i), sweepI_j(i)) = v_i_t(4); 
                  utemp(5, sweepI_i(i), sweepI_j(i)) = v_i_t(5); 
                  utemp(1, sweepI_i(i - 1), sweepI_j(i - 1)) = v_p_t(1); 
                  utemp(2, sweepI_i(i - 1), sweepI_j(i - 1)) = v_p_t(2); 
                  utemp(3, sweepI_i(i - 1), sweepI_j(i - 1)) = v_p_t(3); 
                  utemp(4, sweepI_i(i - 1), sweepI_j(i - 1)) = v_p_t(4); 
                  utemp(5, sweepI_i(i - 1), sweepI_j(i - 1)) = v_p_t(5); 

                  pres_v_p = gamma1*(v_p(4) - 0.5*v_p(2)*v_p(2)/v_p(1) - 0.5*v_p(3)**2/v_p(1) - q_const*v_p(5))
                  pres_v_i_t = gamma1*(v_i_t(4) - 0.5*v_i_t(2)**2/v_i_t(1) - 0.5*v_i_t(3)**2/v_i_t(1) - q_const*v_i_t(5))
                  pres_v_p_t = gamma1*(v_p_t(4) - 0.5*v_p_t(2)**2/v_p_t(1) - 0.5*v_p_t(3)**2/v_p_t(1) - q_const*v_p_t(5))

                  if (pres_v_i + pres_v_p > pres_v_i_t + pres_v_p_t + eps) then
                     ! Distance between v_p and v_i
                     dist_pi = sqrt(sum((v_p - v_i)**2))

                     ! Distance between v_p_t and v_i_t
                     dist_tilde = sqrt(sum((v_p_t - v_i_t)**2))

                     stop "Sweep check failed: total pressure increased after adjustment"
                     ! stop
                  end if

               end if
            end do

         else

            !Forward Sweep
            DO i = 1, (total_x*total_y - 1)

               v_i(1) = utemp(1, sweepII_i(i), sweepII_j(i)); 
               v_i(2) = utemp(2, sweepII_i(i), sweepII_j(i)); 
               v_i(3) = utemp(3, sweepII_i(i), sweepII_j(i)); 
               v_i(4) = utemp(4, sweepII_i(i), sweepII_j(i)); 
               v_i(5) = utemp(5, sweepI_i(i), sweepI_j(i)); 
               pres_v_i = gamma1*(v_i(4) - 0.5*v_i(2)*v_i(2)/v_i(1) - 0.5*v_i(3)**2/v_i(1) - q_const*v_i(5))

               w_i = weights(sweepII_i(i), sweepII_j(i)); 
               !if (pres_v_i <= 0.0 .OR. v_i(1) <= 0.0) then
               if (pres_v_i <= eps .OR. v_i(1) <= eps) then

                  v_p(1) = utemp(1, sweepII_i(i + 1), sweepII_j(i + 1)); 
                  v_p(2) = utemp(2, sweepII_i(i + 1), sweepII_j(i + 1)); 
                  v_p(3) = utemp(3, sweepII_i(i + 1), sweepII_j(i + 1)); 
                  v_p(4) = utemp(4, sweepII_i(i + 1), sweepII_j(i + 1)); 
                  v_p(5) = utemp(5, sweepI_i(i - 1), sweepI_j(i - 1)); 
                  w_p = weights(sweepII_i(i + 1), sweepII_j(i + 1)); 
                  v_i_t = nod_lim(v_i, v_p, u_bar); 
                  d = v_i - v_i_t; 
                  v_p_t = v_p + (w_i/w_p)*(d); 
                  utemp(1, sweepII_i(i), sweepII_j(i)) = v_i_t(1); 
                  utemp(2, sweepII_i(i), sweepII_j(i)) = v_i_t(2); 
                  utemp(3, sweepII_i(i), sweepII_j(i)) = v_i_t(3); 
                  utemp(4, sweepII_i(i), sweepII_j(i)) = v_i_t(4); 
                  utemp(5, sweepI_i(i), sweepI_j(i)) = v_i_t(5); 
                  utemp(1, sweepII_i(i + 1), sweepII_j(i + 1)) = v_p_t(1); 
                  utemp(2, sweepII_i(i + 1), sweepII_j(i + 1)) = v_p_t(2); 
                  utemp(3, sweepII_i(i + 1), sweepII_j(i + 1)) = v_p_t(3); 
                  utemp(4, sweepII_i(i + 1), sweepII_j(i + 1)) = v_p_t(4); 
                  utemp(5, sweepI_i(i - 1), sweepI_j(i - 1)) = v_p_t(5); 

                  pres_v_p = gamma1*(v_p(4) - 0.5*v_p(2)*v_p(2)/v_p(1) - 0.5*v_p(3)**2/v_p(1) - q_const*v_p(5))
                  pres_v_i_t = gamma1*(v_i_t(4) - 0.5*v_i_t(2)**2/v_i_t(1) - 0.5*v_i_t(3)**2/v_i_t(1) - q_const*v_i_t(5))
                  pres_v_p_t = gamma1*(v_p_t(4) - 0.5*v_p_t(2)**2/v_p_t(1) - 0.5*v_p_t(3)**2/v_p_t(1) - q_const*v_p_t(5))

                  if (pres_v_i + pres_v_p > pres_v_i_t + pres_v_p_t + eps) then

                     ! Distance between v_p and v_i
                     dist_pi = sqrt(sum((v_p - v_i)**2))

                     ! Distance between v_p_t and v_i_t
                     dist_tilde = sqrt(sum((v_p_t - v_i_t)**2))

                     !  stop
                  end if

                  ! if (i == 66 .AND. y_iter == 64) then
                  !end if
               end if
            end do

            !Backward Sweep
            DO i = total_x*total_y, 2, -1

               v_i(1) = utemp(1, sweepII_i(i), sweepII_j(i))
               v_i(2) = utemp(2, sweepII_i(i), sweepII_j(i))
               v_i(3) = utemp(3, sweepII_i(i), sweepII_j(i))
               v_i(4) = utemp(4, sweepII_i(i), sweepII_j(i))
               v_i(5) = utemp(5, sweepI_i(i), sweepI_j(i)); 
               pres_v_i = gamma1*(v_i(4) - 0.5*v_i(2)*v_i(2)/v_i(1) - 0.5*v_i(3)**2/v_i(1) - q_const*v_i(5))

               w_i = weights(sweepII_i(i), sweepII_j(i)); 
               if (pres_v_i <= eps .OR. v_i(1) <= eps) then

                  v_p(1) = utemp(1, sweepII_i(i - 1), sweepII_j(i - 1)); 
                  v_p(2) = utemp(2, sweepII_i(i - 1), sweepII_j(i - 1)); 
                  v_p(3) = utemp(3, sweepII_i(i - 1), sweepII_j(i - 1)); 
                  v_p(4) = utemp(4, sweepII_i(i - 1), sweepII_j(i - 1)); 
                  v_p(5) = utemp(5, sweepI_i(i - 1), sweepI_j(i - 1)); 
                  w_p = weights(sweepII_i(i - 1), sweepII_j(i - 1)); 
                  v_i_t = nod_lim(v_i, v_p, u_bar); 
                  d = v_i - v_i_t; 
                  v_p_t = v_p + (w_i/w_p)*(d); 
                  utemp(1, sweepII_i(i), sweepII_j(i)) = v_i_t(1); 
                  utemp(2, sweepII_i(i), sweepII_j(i)) = v_i_t(2); 
                  utemp(3, sweepII_i(i), sweepII_j(i)) = v_i_t(3); 
                  utemp(4, sweepII_i(i), sweepII_j(i)) = v_i_t(4); 
                  utemp(5, sweepI_i(i), sweepI_j(i)) = v_i_t(5); 
                  utemp(1, sweepII_i(i - 1), sweepII_j(i - 1)) = v_p_t(1); 
                  utemp(2, sweepII_i(i - 1), sweepII_j(i - 1)) = v_p_t(2); 
                  utemp(3, sweepII_i(i - 1), sweepII_j(i - 1)) = v_p_t(3); 
                  utemp(4, sweepII_i(i - 1), sweepII_j(i - 1)) = v_p_t(4); 
                  utemp(5, sweepI_i(i - 1), sweepI_j(i - 1)) = v_p_t(5); 

                  pres_v_p = gamma1*(v_p(4) - 0.5*v_p(2)*v_p(2)/v_p(1) - 0.5*v_p(3)**2/v_p(1) - q_const*v_p(5))
                  pres_v_i_t = gamma1*(v_i_t(4) - 0.5*v_i_t(2)**2/v_i_t(1) - 0.5*v_i_t(3)**2/v_i_t(1) - q_const*v_i_t(5))
                  pres_v_p_t = gamma1*(v_p_t(4) - 0.5*v_p_t(2)**2/v_p_t(1) - 0.5*v_p_t(3)**2/v_p_t(1) - q_const*v_p_t(5))

                  if (pres_v_i + pres_v_p > pres_v_i_t + pres_v_p_t + eps) then

                     ! Distance between v_p and v_i
                     dist_pi = sqrt(sum((v_p - v_i)**2))

                     ! Distance between v_p_t and v_i_t
                     dist_tilde = sqrt(sum((v_p_t - v_i_t)**2))

                     !  stop
                  end if

               end if
            end do
         end if

         rho_min = 1d20
         P_min = 1d20
         DO j = 1, total_y
            DO i = 1, total_x
               !P_check(j) = gamma1*(utemp(4,x_iter,j) - 0.5*(utemp(2,x_iter,j)**2 + utemp(3,x_iter,j)**2)/utemp(1,x_iter,j));
    P_min = min(gamma1*(utemp(4, i, j) - 0.5*(utemp(2, i, j)*utemp(2, i, j) + utemp(3, i, j)*utemp(3, i, j))/utemp(1, i, j)), P_min); 
               rho_min = min(utemp(1, i, j), rho_min); 
            end do
         end do

         pre_count = pre_count + 1

         if (P_min > eps .AND. rho_min > eps) then !(P_min > 0 .AND. rho_min>0) then
            sweep_track(sweep_count) = q
            neg_track(sweep_count) = t_neg_count
            pos_track(sweep_count) = t_pos_count

            exit
         end if

         if (q == swp_num) then
            print *, "SWEEP FAIL: swp_num =", swp_num
            stop
         end if

      end do
      !end if

      postag = 0

   END SUBROUTINE sweep_lim

!*****************************
   SUBROUTINE scalar_sweep
      IMPLICIT NONE
      REAL(8)::  v_j, v_jp1, v_jm1, v_i, v_ip1, v_im1
      REAL(8), dimension(1:total_x, 1:total_y):: rho_t
      REAL(8)::  w_j, w_jp1, w_jm1, w_i, w_ip1, w_im1, delta, rho_min, rho_max
      INTEGER:: ss

      if (use_limiter == 2) then
         rho_min = 1d20
         DO j = 1, total_y
            DO i = 1, total_x
               rho_min = min(utemp(1, i, j), rho_min); 
            end do
         end do

         delta = 1d-12; 
         if (rho_min < delta) then
         DO j = 1, total_y
            DO i = 1, total_x
               rho_t(i, j) = utemp(1, i, j); 
            end do
         end do

         DO i = 1, (total_x*total_y - 1)
            v_i = rho_t(sweepI_i(i), sweepI_j(i)); v_ip1 = rho_t(sweepI_i(i + 1), sweepI_j(i + 1)); w_i = weights(sweepI_i(i), sweepI_j(i)); w_ip1 = weights(sweepI_i(i + 1), sweepI_j(i + 1)); 
            if (v_i < delta) then
               rho_t(sweepI_i(i + 1), sweepI_j(i + 1)) = v_ip1 + (w_i/w_ip1)*(v_i - delta); 
               rho_t(sweepI_i(i), sweepI_j(i)) = delta; 
            end if
         end do

         DO i = total_x*total_y, 2, -1
            v_i = rho_t(sweepI_i(i), sweepI_j(i)); v_im1 = rho_t(sweepI_i(i - 1), sweepI_j(i - 1)); w_i = weights(sweepI_i(i), sweepI_j(i)); w_im1 = weights(sweepI_i(i - 1), sweepI_j(i - 1)); 
            if (v_i < delta) then
               rho_t(sweepI_i(i - 1), sweepI_j(i - 1)) = v_im1 + (w_i/w_im1)*(v_i - delta); 
               rho_t(sweepI_i(i), sweepI_j(i)) = delta; 
            end if
         end do

         !Fill in values for u!
         DO j = 1, total_y
            DO i = 1, total_x
               utemp(1, i, j) = rho_t(i, j); 
            end do
         end do

         den_count = den_count + 1; 
         end if
         rho_min = 1d20
         DO j = 1, total_y
            DO i = 1, total_x
               rho_min = min(utemp(1, i, j), rho_min); 
            end do
         end do

         rho_max = 0
         DO j = 1, total_y
            DO i = 1, total_x
               rho_max = max(utemp(1, i, j), rho_max); 
            end do
         end do

         postag = 0

      end if

   END SUBROUTINE scalar_sweep

   SUBROUTINE sweeping_order
      IMPLICIT NONE

! ====================
! Sweep I: column-by-column
! ====================
      k = 0
      do j = 1, total_y
      if (mod(j, 2) == 1) then
         ! odd column: sweep upward (i=1..total_x)
         do i = 1, total_x
            k = k + 1
            sweepI_i(k) = i
            sweepI_j(k) = j
         end do
      else
! even column: sweep downward (i=total_x..1)
         do i = total_x, 1, -1
            k = k + 1
            sweepI_i(k) = i
            sweepI_j(k) = j
         end do
      end if
      end do

! ====================
! Sweep II: row-by-row
! ====================
      k = 0
      do i = 1, total_x
      if (mod(i, 2) == 1) then
! odd row: sweep left-to-right (j=1..total_y)
         do j = 1, total_y
            k = k + 1
            sweepII_i(k) = i
            sweepII_j(k) = j
         end do
      else
! even row: sweep right-to-left (j=total_y..1)
         do j = total_y, 1, -1
            k = k + 1
            sweepII_i(k) = i
            sweepII_j(k) = j
         end do
      end if
      end do

!do k = 1, total_x*total_y
!end do

!do k = 1, total_x*total_y
!end do

!deallocate(sweepI_i, sweepI_j, sweepII_i, sweepII_j)
   END SUBROUTINE sweeping_order

   subroutine gather_to_root
      implicit none
      REAL(8):: recv_buf(5, 1:nx, 1:ny)
      REAL(8), ALLOCATABLE :: buffer(:)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: buffer_addr
      INTEGER :: buffer_size_bytes
      !INTEGER :: mpierror

      real(8):: lrho, lm, ln, le, pres
      integer:: den_tmp, pre_tmp

! Calculate the size of the buffer in bytes
      ! (size of REAL(8) is 8 bytes)
      buffer_size_bytes = (5*(cut_x*nx + 6)*(cut_y*ny + 6))*SIZEOF(1.0_8)

! Allocate the buffer
      ALLOCATE (buffer(buffer_size_bytes))

      ! Assumes these variables are already declared in the program:
      ! u, utemp, nx, ny, cut_x, cut_y, mpirank, mpisize, recv_buf, send_status
      ! buffer, buffer_size_bytes, mpierror, mpitag, mpides, mpisource

      ! Attach buffer for buffered sends
      CALL MPI_BUFFER_ATTACH(buffer, buffer_size_bytes, mpierror)

      if (mpirank .ne. 0) then
         mpitag = mpirank
         mpides = 0
         CALL MPI_SEND(u(:, 1:nx, 1:ny), 5*nx*ny, MPI_REAL8, mpides, mpitag, MPI_COMM_WORLD, mpierror)
      else
         ! Root copies its own block
         DO ig = 1, nx
            DO jg = 1, ny
               utemp(:, ig, jg) = u(:, ig, jg)
            END DO
         END DO

         ! Receive from all other ranks
         DO i = 1, mpisize - 1
            mpisource = i
            mpitag = i
            CALL MPI_RECV(recv_buf, 5*nx*ny, MPI_REAL8, mpisource, mpitag, MPI_COMM_WORLD, send_status, mpierror)

            DO ig = 1, nx
               DO jg = 1, ny
                  utemp(:, ig + mod(i, cut_x)*nx, jg + i/cut_x*ny) = recv_buf(:, ig, jg)
               END DO
            END DO
         END DO
      END IF

      ! Detach and deallocate buffer
      CALL MPI_BUFFER_DETACH(buffer, buffer_size_bytes, mpierror)
      DEALLOCATE (buffer)
   end subroutine gather_to_root

!-------------------------------------------------------------
   subroutine scatter_from_root
      implicit none
      REAL(8):: recv_buf(5, 1:nx, 1:ny)
      REAL(8) :: send_buf(5, 1:nx, 1:ny)

      REAL(8), ALLOCATABLE :: buffer(:)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: buffer_addr
      INTEGER :: buffer_size_bytes
      !INTEGER :: mpierror

      real(8):: lrho, lm, ln, le, pres
      integer:: den_tmp, pre_tmp

! Calculate the size of the buffer in bytes
      ! (size of REAL(8) is 8 bytes)
      buffer_size_bytes = (5*(cut_x*nx + 6)*(cut_y*ny + 6))*SIZEOF(1.0_8)

! Allocate the buffer
      ALLOCATE (buffer(buffer_size_bytes))

      ! Assumes the same variables as gather:
      ! u, utemp, nx, ny, cut_x, cut_y, mpirank, mpisize, recv_buf, send_buf, recv_status
      ! buffer, buffer_size_bytes, mpierror, mpitag, mpides, mpisource

      ! Attach buffer for buffered sends
      CALL MPI_BUFFER_ATTACH(buffer, buffer_size_bytes, mpierror)

      if (mpirank == 0) then
         ! Root copies its own block
         DO ig = 1, nx
            DO jg = 1, ny
               u(:, ig, jg) = utemp(:, ig, jg)
            END DO
         END DO

         ! Send to all other ranks
         DO i = 1, mpisize - 1
            mpides = i
            mpitag = i

            DO ig = 1, nx
               DO jg = 1, ny
                  send_buf(:, ig, jg) = utemp(:, ig + mod(i, cut_x)*nx, jg + i/cut_x*ny)
               END DO
            END DO

            CALL MPI_BSEND(send_buf, 5*nx*ny, MPI_REAL8, mpides, mpitag, MPI_COMM_WORLD, mpierror)
         END DO

      else
         ! Non-root receives its block
         mpisource = 0
         mpitag = mpirank
         CALL MPI_RECV(recv_buf, 5*nx*ny, MPI_REAL8, mpisource, mpitag, MPI_COMM_WORLD, recv_status, mpierror)

         DO ig = 1, nx
            DO jg = 1, ny
               u(:, ig, jg) = recv_buf(:, ig, jg)
            END DO
         END DO
      END IF

      ! Detach and deallocate buffer
      CALL MPI_BUFFER_DETACH(buffer, buffer_size_bytes, mpierror)
      DEALLOCATE (buffer)
   end subroutine scatter_from_root
END PROGRAM FDWENO4EULER1D
