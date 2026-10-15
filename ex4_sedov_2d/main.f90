PROGRAM FDWENO4EULER1D
   USE weno5
   use MPI
   IMPLICIT NONE

!********************************************************
! This is a parallel MPI code for 3rd order RK + 5th order
! finite difference WENO for the 2D Euler equation.
!                     Xiangxiong Zhang  Brown University
!                                         July 2011
!  Edited to test sweeping algorithm
!                      Chloe Griffin Brown University
!                                         October 2025
!********************************************************
   INTEGER, PARAMETER::  ierror = 0, mg = 0
   
   INTEGER, PARAMETER :: cut_x = 8, cut_y = 4
   INTEGER, PARAMETER :: total_x = 640, total_y =  640
   INTEGER, PARAMETER :: nx = total_x/cut_x, ny = total_y / cut_y
   REAL(8), PARAMETER:: time_terminal=1, gamma=1.4, gamma1=gamma-1.,  eps=1d-13, eps2=1d-8, xm=1.3, ym=1.3
   REAL(8), DIMENSION(4, -2:nx + 3, -2:ny + 3):: u, u0, RHS, fp_flux, fn_flux, fp, fn, gp_flux, gn_flux, gp, gn
   REAL(8), DIMENSION(-2:nx + 3):: x
   REAL(8), DIMENSION(-2:ny + 3):: y
   REAL(8), DIMENSION(-2:nx + 3, -2:ny + 3)::alpha
   integer:: rktag, kcount
   real(8):: dx, dy, dt, time, alphamax, xlen, ylen
   integer:: ig, jg, i, j, k, kk, m, den_count = 0, pre_count = 0, postag
   
   !REAL(8):: utemp(4, 1:total_x, 1:total_y)
   
   !REAL(8), DIMENSION(4, 1:total_x, 1:total_y):: utemp = 1.0
 
   real(8), allocatable:: utemp(:,:,:)


   INTEGER:: mpirank, mpierror, mpisize, position_x, position_y, destination(8), source(8)&
  &, mpitag, mpisource, recv_status(MPI_STATUS_SIZE), send_status(MPI_STATUS_SIZE), mpides&
  &, send_request(4), recv_request(4), color_x, color_y
   REAL:: time_start, time_end
   real:: xm1, ss, dri, uri, vri, pri, ach2, ple, ule, dle, vle

   !Adding this for dynamics output storage
   character(len=256) :: output_dir, err_file_path, den_file_path, m1_file_path
   character(len=256) :: m2_file_path, pre_file_path
   integer :: ios

   !Adding sweeping types
    !Adding sweeping types
   integer::  x_iter, y_iter, iter
   INTEGER, PARAMETER :: use_limiter = 2, swp_num = 20000
   Real(8), PARAMETER:: nSteps_x = total_x, nsteps_y = total_y

   INTEGER, PARAMETER :: local_size = nx * ny
   REAL(8), DIMENSION(1:total_x, 1:total_y):: weights
   integer, allocatable :: sweepI_i(:), sweepI_j(:)
   integer, allocatable :: sweepII_i(:), sweepII_j(:)

    !Adding sweep data
   INTEGER, PARAMETER :: total_max = 5000
   INTEGER, DIMENSION(1:total_max) :: sweep_track = 0
   INTEGER :: sweep_count = 0, total_sweep, max_sweep
   real(8):: avg_sweep

   INTEGER, DIMENSION(1:total_max) :: neg_track = 0
   INTEGER, DIMENSION(1:total_max) :: pos_track = 0
   INTEGER :: t_neg_count, t_pos_count, total_neg, total_pos, max_neg, max_pos
   real(8):: avg_pos, avg_neg



   call MPI_Init(mpierror)
   call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpierror)
   call MPI_Comm_size(MPI_COMM_WORLD, mpisize, mpierror)

    if (mpirank == 0) then
      allocate(utemp(4, 1:total_x, 1:total_y))
   end if 

   allocate(sweepI_i(total_x*total_y), sweepI_j(total_x*total_y))
   allocate(sweepII_i(total_x*total_y), sweepII_j(total_x*total_y))

   if (mpirank == 0) then
      !allocate(sweepI_i(total_x*total_y), sweepI_j(total_x*total_y))
      !allocate(sweepII_i(total_x*total_y), sweepII_j(total_x*total_y))
      call sweeping_order
   end if

   ! Get command-line arguments
   call get_command_argument(1, output_dir)

   
   IF (mpisize .ne. cut_x*cut_y) then
      write (*, *) 'Need ', cut_x*cut_y, ' nodes!!!'
      stop
   END IF

   ! Set the local grid size for each process
   !nx = total_x / cut_x ! This will be 400 / 8 = 50
   !ny = total_y / cut_y ! This will be 800 / 4 = 200

     DO j = 1, total_y
      DO i = 1, total_x 
      weights(i,j) = 1./(nSteps_x*nSteps_y)
      END DO
   END DO

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
      REAL(8), INTENT(IN):: a(4)
      REAL(8):: vx, vy, pres, rho, mon_x, mon_y, e
      REAL(8), DIMENSION(1:4)::fflux

      rho = a(1)
      mon_x = a(2)
      mon_y = a(3)
      e = a(4)

      vx = mon_x/rho
      vy = mon_y/rho

      if (rho < eps2) vx = 0.
      if (rho < eps2) vy = 0.

      pres = gamma1*(e - 0.5*rho*vx**2 - 0.5*rho*vy**2)

      fflux(1) = mon_x
      fflux(2) = rho*vx**2 + pres
      fflux(3) = rho*vx*vy
      fflux(4) = (e + pres)*vx

      RETURN
   END FUNCTION fflux

   FUNCTION gflux(a)
      IMPLICIT NONE
      REAL(8), INTENT(IN):: a(4)
      REAL(8):: vx, vy, pres, rho, mon_x, mon_y, e
      REAL(8), DIMENSION(1:4)::gflux

      rho = a(1)
      mon_x = a(2)
      mon_y = a(3)
      e = a(4)

      vx = mon_x/rho
      vy = mon_y/rho

      if (rho < eps2) vx = 0.
      if (rho < eps2) vy = 0.

      pres = gamma1*(e - 0.5*rho*vx**2 - 0.5*rho*vy**2)

      gflux(1) = mon_y
      gflux(2) = rho*vx*vy
      gflux(3) = rho*vy**2 + pres
      gflux(4) = (e + pres)*vy

      RETURN
   END FUNCTION gflux
   !**************************
   SUBROUTINE wenofd
      IMPLICIT NONE
      REAL(8):: evr(4, 4), evl(4, 4), cu_p(4, -2:2), cu_n(4, -2:2)
      REAL(8), DIMENSION(0:max(nx, ny) + 2):: w, vx, vy, h
      REAL(8):: den, xmt, ymt, eng, t0, t1, t2, t3, pre, vex, vey, cvel, vxm, vym, hm, qm, cm, rcm, b1, b2
      REAL(8):: ul(4), ur(4), d_ave, m_ave, e_ave, p_ave, theta, point(4, 2), d_l, m_l, e_l, p_l, n_l, n_ave
      real(8):: alpha_local

      call findalpha

      !do i=-2,nx+3
      !do j=-2,ny+3

      !  fp(:,i,j)=0.5*(u(:,i,j)+fflux(u(:,i,j))/alphamax)
      !fn(:,i,j)=0.5*(u(:,i,j)-fflux(u(:,i,j))/alphamax)
      !gp(:,i,j)=0.5*(u(:,i,j)+gflux(u(:,i,j))/alphamax)
      !gn(:,i,j)=0.5*(u(:,i,j)-gflux(u(:,i,j))/alphamax)
      !enddo
      !enddo

      DO j = 0, ny + 1

         do i = 0, nx + 2
            den = u(1, i, j)
            xmt = u(2, i, j)
            ymt = u(3, i, j)
            eng = u(4, i, j)
            t0 = 1./den
            vex = xmt*t0
            vey = ymt*t0
            pre = gamma1*(eng - 0.5*(xmt*vex + ymt*vey))
            w(i) = sqrt(abs(den))
            vx(i) = vex
            vy(i) = vey
            h(i) = (pre + eng)*t0
         end do

         alpha_local = 0.0_8
         do i = -2, nx + 3
            do ig = -mg, mg
               alpha_local = max(alpha_local, alpha(i, j + ig))
            end do
         end do

         alpha_local = alphamax
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
            hm = t0*h(i) + t1*h(i + 1)
            qm = 0.5*(vxm*vxm + vym*vym)
            cm = sqrt(abs(gamma1*(hm - qm))) + 1d-15
            t0 = vxm*cm
            evr(1, 1) = 1.0
            evr(1, 2) = 0.0
            evr(1, 3) = 1.0
            evr(1, 4) = 1.0
            evr(2, 1) = vxm - cm
            evr(2, 2) = 0.0
            evr(2, 3) = vxm
            evr(2, 4) = vxm + cm
            evr(3, 1) = vym
            evr(3, 2) = 1.0
            evr(3, 3) = vym
            evr(3, 4) = vym
            evr(4, 1) = hm - t0
            evr(4, 2) = vym
            evr(4, 3) = qm
            evr(4, 4) = hm + t0
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
            evl(2, 1) = -vym
            evl(2, 2) = 0.0
            evl(2, 3) = 1.0
            evl(2, 4) = 0.0
            evl(3, 1) = 1.-b2
            evl(3, 2) = t1
            evl(3, 3) = t3
            evl(3, 4) = -b1
            evl(4, 1) = 0.5*(b2 - t0)
            evl(4, 2) = -0.5*(t1 - rcm)
            evl(4, 3) = -0.5*t3
            evl(4, 4) = t2

            !******************************************
            ! projection into characteristic fields
            !******************************************
            DO ig = -2, 2
            DO k = 1, 4
     cu_p(k, ig) = evl(k, 1)*fp(1, i + ig, j) + evl(k, 2)*fp(2, i + ig, j) + evl(k, 3)*fp(3, i + ig, j) + evl(k, 4)*fp(4, i + ig, j)
     cu_n(k, ig) = evl(k, 1)*fn(1, i + ig, j) + evl(k, 2)*fn(2, i + ig, j) + evl(k, 3)*fn(3, i + ig, j) + evl(k, 4)*fn(4, i + ig, j)
            END DO
            END DO
            !**************************
            ! weno reconstruction
            !**************************
            DO kk = 1, 4
               call weno_right(cu_p(kk, :), ur(kk))
               call weno_left(cu_n(kk, :), ul(kk))
            END DO
            !******************************************
            ! projection back into conserved variables
            !******************************************
            DO k = 1, 4
               fp_flux(k, i, j) = evr(k, 1)*ur(1) + evr(k, 2)*ur(2) + evr(k, 3)*ur(3) + evr(k, 4)*ur(4)
               fn_flux(k, i, j) = evr(k, 1)*ul(1) + evr(k, 2)*ul(2) + evr(k, 3)*ul(3) + evr(k, 4)*ul(4)
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
            p_ave = gamma1*(e_ave - 0.5*m_ave**2/d_ave - 0.5*n_ave**2/d_ave)
            theta = 1.
            point(:, 1) = fp_flux(:, i, j)
            point(:, 2) = (fp(:, i, j) - fp_flux(:, i, j)/12.)/(11./12.)

            do ig = 1, 2
               d_l = point(1, ig)
               if (d_l < 0.) theta = min(theta, abs(d_ave - 0.)/abs(d_ave - d_l))
            end do
           ! theta = max(theta, 0.) !This is not in other code

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

               p_l = gamma1*(e_l - m_l**2/d_l/2.-n_l**2/d_l/2.)
               if (p_l < 0.) theta = min(theta, abs(p_ave - 0.)/abs(p_ave - p_l))
            end do

            !theta = max(theta, 0.)
            ! if (u(1,i,j)<1d-3) theta=0.
            if (theta < 1d0) then
               pre_count = pre_count + 1
               fp_flux(1, i, j) = theta*fp_flux(1, i, j) + (1.-theta)*d_ave
               fp_flux(2, i, j) = theta*fp_flux(2, i, j) + (1.-theta)*m_ave
               fp_flux(3, i, j) = theta*fp_flux(3, i, j) + (1.-theta)*n_ave
               fp_flux(4, i, j) = theta*fp_flux(4, i, j) + (1.-theta)*e_ave

            else
            end if

            theta = 1d0
            d_ave = fn(1, i, j)
            m_ave = fn(2, i, j)
            n_ave = fn(3, i, j)
            e_ave = fn(4, i, j)
            p_ave = gamma1*(e_ave - 0.5*m_ave**2/d_ave - 0.5*n_ave**2/d_ave)
            theta = 1.
            point(:, 1) = fn_flux(:, i, j)
            point(:, 2) = (fn(:, i, j) - fn_flux(:, i, j)/12.)/(11./12.)

            do ig = 1, 2
               d_l = point(1, ig)
               if (d_l < 0.) theta = min(theta, abs(d_ave - 0.)/abs(d_ave - d_l))
            end do

           ! theta = max(theta, 0.)
            ! if (u(1,i,j)<1d-3) theta=0.
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
               p_l = gamma1*(e_l - m_l**2/d_l/2.-n_l**2/d_l/2.)
               if (p_l < 0.) theta = min(theta, abs(p_ave - 0.)/abs(p_ave - p_l))
            end do
            !theta = max(theta, 0.)
            ! if (u(1,i,j)<1d-3) theta=0.
            if (theta < 1d0) then
               pre_count = pre_count + 1
               fn_flux(1, i, j) = theta*fn_flux(1, i, j) + (1.-theta)*d_ave
               fn_flux(2, i, j) = theta*fn_flux(2, i, j) + (1.-theta)*m_ave
               fn_flux(3, i, j) = theta*fn_flux(3, i, j) + (1.-theta)*n_ave
               fn_flux(4, i, j) = theta*fn_flux(4, i, j) + (1.-theta)*e_ave
            else
            end if
         end if

            fp_flux(:, i, j) = alpha_local*fp_flux(:, i, j) !Other code multiply later
            fn_flux(:, i, j) = alpha_local*fn_flux(:, i, j)
         END DO
      END DO

      DO i = 0, nx + 1

         do j = 0, ny + 2
            den = u(1, i, j)
            xmt = u(2, i, j)
            ymt = u(3, i, j)
            eng = u(4, i, j)
            t0 = 1./(den + 1d-15)
            vex = xmt*t0
            vey = ymt*t0
            pre = gamma1*(eng - 0.5*(xmt*vex + ymt*vey))
            w(j) = sqrt(abs(den))
            vx(j) = vex
            vy(j) = vey
            h(j) = (pre + eng)*t0
         end do

         alpha_local = 0.0_8 !Different than other 2D Sedov code
         do j = -2, ny + 3
            do ig = -mg, mg
               alpha_local = max(alpha_local, alpha(i + ig, j))
            end do
         end do
         !alpha_local=5*alpha_local
         alpha_local = alphamax
         do j = -2, ny + 3 !Here doing what on top in other code
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
            qm = 0.5*(vxm*vxm + vym*vym)
            cm = sqrt(abs(gamma1*(hm - qm))) + 1d-15
            t0 = vym*cm
            evr(1, 1) = 1.0
            evr(1, 2) = 0.0
            evr(1, 3) = 1.0
            evr(1, 4) = 1.0
            evr(2, 1) = vxm
            evr(2, 2) = 1.0
            evr(2, 3) = vxm
            evr(2, 4) = vxm
            evr(3, 1) = vym - cm
            evr(3, 2) = 0.0
            evr(3, 3) = vym
            evr(3, 4) = vym + cm
            evr(4, 1) = hm - t0
            evr(4, 2) = vxm
            evr(4, 3) = qm
            evr(4, 4) = hm + t0
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
            evl(2, 1) = -vxm
            evl(2, 2) = 1.0
            evl(2, 3) = 0.0
            evl(2, 4) = 0.0
            evl(3, 1) = 1.-b2
            evl(3, 2) = t3
            evl(3, 3) = t1
            evl(3, 4) = -b1
            evl(4, 1) = 0.5*(b2 - t0)
            evl(4, 2) = -0.5*t3
            evl(4, 3) = -0.5*(t1 - rcm)
            evl(4, 4) = t2

            !******************************************
            ! projection into characteristic fields
            !******************************************
            DO ig = -2, 2
            DO k = 1, 4
     cu_p(k, ig) = evl(k, 1)*gp(1, i, j + ig) + evl(k, 2)*gp(2, i, j + ig) + evl(k, 3)*gp(3, i, j + ig) + evl(k, 4)*gp(4, i, j + ig)
     cu_n(k, ig) = evl(k, 1)*gn(1, i, j + ig) + evl(k, 2)*gn(2, i, j + ig) + evl(k, 3)*gn(3, i, j + ig) + evl(k, 4)*gn(4, i, j + ig)
            END DO
            END DO
            !**************************
            ! weno reconstruction
            !**************************
            DO kk = 1, 4
               call weno_right(cu_p(kk, :), ur(kk))
               call weno_left(cu_n(kk, :), ul(kk))
            END DO
            !******************************************
            ! projection back into conserved variables
            !******************************************
            DO k = 1, 4
               gp_flux(k, i, j) = evr(k, 1)*ur(1) + evr(k, 2)*ur(2) + evr(k, 3)*ur(3) + evr(k, 4)*ur(4)
               gn_flux(k, i, j) = evr(k, 1)*ul(1) + evr(k, 2)*ul(2) + evr(k, 3)*ul(3) + evr(k, 4)*ul(4)
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
            p_ave = gamma1*(e_ave - 0.5*m_ave**2/d_ave - 0.5*n_ave**2/d_ave)
            theta = 1.
            point(:, 1) = gp_flux(:, i, j)
            point(:, 2) = (gp(:, i, j) - gp_flux(:, i, j)/12.)/(11./12.)

            do ig = 1, 2
               d_l = point(1, ig)
               if (d_l < 0.) theta = min(theta, abs(d_ave - 0.)/abs(d_ave - d_l))
            end do
            !theta = max(theta, 0.)
            ! if (u(1,i,j)<1d-3) theta=0.
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
               p_l = gamma1*(e_l - m_l**2/d_l/2.-n_l**2/d_l/2.)
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

            else
            end if

            theta = 1d0
            d_ave = gn(1, i, j)
            m_ave = gn(2, i, j)
            n_ave = gn(3, i, j)
            e_ave = gn(4, i, j)
            p_ave = gamma1*(e_ave - 0.5*m_ave**2/d_ave - 0.5*n_ave**2/d_ave)
            theta = 1.
            point(:, 1) = gn_flux(:, i, j)
            point(:, 2) = (gn(:, i, j) - gn_flux(:, i, j)/12.)/(11./12.)

            do ig = 1, 2
               d_l = point(1, ig)
               if (d_l < 0.) theta = min(theta, abs(d_ave - 0.)/abs(d_ave - d_l))
            end do
            !theta = max(theta, 0.)
            ! if (u(1,i,j)<1d-3) theta=0.
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
               p_l = gamma1*(e_l - m_l**2/d_l/2.-n_l**2/d_l/2.)
               if (p_l < 0.) theta = min(theta, abs(p_ave - 0.)/abs(p_ave - p_l))
            end do
            !theta = max(theta, 0.)
            ! if (u(1,i,j)<1d-3) theta=0.
            if (theta < 1d0) then
               pre_count = pre_count + 1
               gn_flux(1, i, j) = theta*gn_flux(1, i, j) + (1.-theta)*d_ave
               gn_flux(2, i, j) = theta*gn_flux(2, i, j) + (1.-theta)*m_ave
               gn_flux(3, i, j) = theta*gn_flux(3, i, j) + (1.-theta)*n_ave
               gn_flux(4, i, j) = theta*gn_flux(4, i, j) + (1.-theta)*e_ave
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
      RHS = 0;
      DO i = 1, nx
      DO j = 1, ny
         RHS(:, i, j) = -((fp_flux(:, i, j) - fn_flux(:, i + 1, j)) - (fp_flux(:, i - 1, j) - fn_flux(:, i, j)))/dx&
         &        - ((gp_flux(:, i, j) - gn_flux(:, i, j + 1)) - (gp_flux(:, i, j - 1) - gn_flux(:, i, j)))/dy

      END DO
      END DO

   END SUBROUTINE scheme
   !**********************
   SUBROUTINE findalpha
      IMPLICIT NONE
      REAL(8)::den, pres, sound, alphatmp1, v_x, v_y
      REAL(8):: alphatmp1_real, alphamax_real 
      
      alphamax=0.
      alpha=0.
      do i=-2,nx+3
      do j=-2,ny+3
         den=u(1,i,j)
      if (den>=eps2) then
         v_x=u(2,i,j)/den
         v_y=u(3,i,j)/den
      else
         v_x=0.
         v_y=0.
      endif
	 
      pres=gamma1*(u(4,i,j)-0.5*den*v_x**2-0.5*den*v_y**2)
      sound=sqrt(gamma*pres/den)

      
      alpha(i,j)=max(abs(u(2,i,j)/u(1,i,j))+sound,abs(u(3,i,j)/u(1,i,j))+sound)
      alphamax=max(alphamax,alpha(i,j))
            
      enddo
      enddo
      call MPI_AllReduce(alphamax, alphatmp1, 1, MPI_REAL8, MPI_MAX, &
      & MPI_COMM_WORLD, mpierror)
      alphamax=alphatmp1
      
   END SUBROUTINE findalpha
   !**************************
   SUBROUTINE rk
      IMPLICIT NONE
      REAL(8), PARAMETER:: cfl = 0.5
      REAL(8):: alphatemp, dt_cfl_x, dt_cfl_y, dt_cfl
      rktag = 1
      CALL wenofd

     ! cfl = 1./2.
     ! dt = dx*cfl/alphamax
      !dt_cfl_x = dx**(5./3.) * cfl / alphamax
      !dt_cfl_y = dy**(5./3.) * cfl / alphamax
      dt_cfl_x = dx * cfl / alphamax
      dt_cfl_y = dy * cfl / alphamax
      dt_cfl = min(dt_cfl_x, dt_cfl_y)

      dt = min(dt_cfl , time_terminal - time)

      alphatemp = alphamax
      postag = 0
      !110 continue
      rktag = 1
      !if (postag==1) then
      !    u=u0
      !dt=dt/2.
      !alphamax=alphatemp
      !else
      !endif

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
         if (u(1, i, j) < eps) then
            postag = 1
         else
         end if

         if (u(4, i, j) - u(3, i, j)**2/u(1, i, j)/2.-u(2, i, j)**2/u(1, i, j)/2. < eps) then
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
     ! integer, parameter:: commsize = 12*(ny + 2 + nx + 2 + 6)*100
     ! real(8), dimension(commsize) :: buffer
!integer, parameter :: commsize = 12*(ny + 2 + nx + 2 + 6)*100
!integer(kind=MPI_ADDRESS_KIND), pointer  :: buffer(:)
!integer(kind=MPI_ADDRESS_KIND) :: buf_size

 INTEGER, PARAMETER :: commsize = 12 * (ny + 2 + nx + 2 + 6) * 100
    
    ! Declare the buffer as a simple real(8) array, NOT a pointer
    REAL(8), DIMENSION(commsize) :: buffer
    
    INTEGER :: buf_size

    REAL(8), DIMENSION(4, 3, -2:ny + 3):: send_buf_y, recv_buf_y
    REAL(8), DIMENSION(4, -2:nx + 3, 3):: send_buf_x, recv_buf_x

integer :: mpierror

!size = commsize

    !  real(8), dimension(4, 3, -2:ny + 3):: send_buf_y, recv_buf_y
    !  real(8), dimension(4, -2:nx + 3, 3):: send_buf_x, recv_buf_x
!buf_size = commsize

buf_size = commsize * 8 ! The size in bytes (real(8) is 8 bytes)

!allocate(buffer(buf_size))
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
         call MPI_BUFFER_ATTACH(buffer, buf_size, mpierror)
         mpitag = 1
         send_buf_y(:, 1, :) = u(:, 1, :)
         send_buf_y(:, 2, :) = u(:, 2, :)
         send_buf_y(:, 3, :) = u(:, 3, :)
         call MPI_Bsend(send_buf_y, 12*(ny + 6), MPI_REAL8,&
        & destination(1), mpitag, MPI_COMM_WORLD, mpierror)
         mpitag = 2
         send_buf_y(:, 1, :) = u(:, nx, :)
         send_buf_y(:, 2, :) = u(:, nx - 1, :)
         send_buf_y(:, 3, :) = u(:, nx - 2, :)
         call MPI_Bsend(send_buf_y, 12*(ny + 6), MPI_REAL8,&
        & destination(2), mpitag, MPI_COMM_WORLD, mpierror)
         mpitag = 3
         send_buf_x(:, :, 1) = u(:, :, ny)
         send_buf_x(:, :, 2) = u(:, :, ny - 1)
         send_buf_x(:, :, 3) = u(:, :, ny - 2)
         call MPI_Bsend(send_buf_x, 12*(nx + 6), MPI_REAL8,&
        & destination(3), mpitag, MPI_COMM_WORLD, mpierror)
         mpitag = 4
         send_buf_x(:, :, 1) = u(:, :, 1)
         send_buf_x(:, :, 2) = u(:, :, 2)
         send_buf_x(:, :, 3) = u(:, :, 3)
         call MPI_Bsend(send_buf_x, 12*(nx + 6), MPI_REAL8,&
        & destination(4), mpitag, MPI_COMM_WORLD, mpierror)

         mpitag = 1
         call MPI_RECV(recv_buf_y, 12*(ny + 6), MPI_REAL8, &
        & source(2), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, nx + 1, :) = recv_buf_y(:, 1, :)
         u(:, nx + 2, :) = recv_buf_y(:, 2, :)
         u(:, nx + 3, :) = recv_buf_y(:, 3, :)
         mpitag = 2
         call MPI_RECV(recv_buf_y, 12*(ny + 6), MPI_REAL8, &
        & source(1), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, 0, :) = recv_buf_y(:, 1, :)
         u(:, -1, :) = recv_buf_y(:, 2, :)
         u(:, -2, :) = recv_buf_y(:, 3, :)

         mpitag = 3
         call MPI_RECV(recv_buf_x, 12*(nx + 6), MPI_REAL8, &
        & source(4), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, :, 0) = recv_buf_x(:, :, 1)
         u(:, :, -1) = recv_buf_x(:, :, 2)
         u(:, :, -2) = recv_buf_x(:, :, 3)
         mpitag = 4
         call MPI_RECV(recv_buf_x, 12*(nx + 6), MPI_REAL8, &
        & source(3), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, :, ny + 1) = recv_buf_x(:, :, 1)
         u(:, :, ny + 2) = recv_buf_x(:, :, 2)
         u(:, :, ny + 3) = recv_buf_x(:, :, 3)

!         call MPI_BUFFER_DETACH(buffer, commsize, mpierror)
         call MPI_BUFFER_DETACH(buffer,buf_size, mpierror) 

      else
      end if

      !************************************************************
      ! outflow for other boundaries
      !************************************************************
      !************************************************************
      ! reflective for x=0, y=0, outflow for other two boundaries
      !************************************************************
      if (position_x == 0) then

         do j = -2, ny + 3
               do i = 1, 3
                    ! i_ghost = 0, -1, -2
               u(1, 1 - i, j) =  u(1, i, j)       ! rho
               u(2, 1 - i, j) = -u(2, i, j)       ! rho*u (normal component reversed)
               u(3, 1 - i, j) =  u(3, i, j)       ! rho*v (tangential component unchanged)
               u(4, 1 - i, j) =  u(4, i, j)       ! energy
               end do
         end do
      else
      end if

      ! Bottom boundary
      if (position_y == 0) then
      do j = 1, 3
      do i = -2, nx+3
          !j_ghost = 0, -1, -2
        u(1, i, 1 - j) =  u(1, i, j)        ! rho
         u(2, i, 1 - j) =  u(2, i, j)       ! rho*u (tangential component unchanged)
         u(3, i, 1 - j) = -u(3, i, j)       ! rho*v (normal component reversed)
         u(4, i, 1 - j) =  u(4, i, j)       ! energy
      end do
      end do
      else
      end if

      ! Right boundary
      !do i = 1, 3
      !do j = -2, ny+3
      !   u(1, nx + i, j) = u(1, nx, j)
      !   u(2, nx + i, j) = u(2, nx, j)
      !   u(3, nx + i, j) = u(3, nx, j)
      !   u(4, nx + i, j) = u(4, nx, j)
      !end do
      !end do

      ! Top boundary (Outflow)
      !do j = 1, 3
      !do i = -2, nx+3
      !   u(1, i, ny + j) = u(1, i, ny)
      !   u(2, i, ny + j) = u(2, i, ny)
      !   u(3, i, ny + j) = u(3, i, ny)
      !   u(4, i, ny + j) = u(4, i, ny)
      !end do
      !end do

      !Old code boundary conditions
      !if (position_y == 0) then
      !   do j = 1, 3
      !   do i = -2, nx + 3
      !     u(:, i, 1 - j) = u(:, i, j)
      !      u(3, i, 1 - j) = -u(3, i, j)
      !   end do
      !   end do
      !else
      !end if


      if (x(nx + 1) > xm) then
         do ig = 1, 3
            u(:, nx + ig, :) = u(:, nx, :)
         end do
      else

      end if
      if (y(ny + 1) > ym) then
         do ig = 1, 3
            u(:, :, ny + ig) = u(:, :, ny)
         end do
      else
      end if

   END SUBROUTINE bc

   SUBROUTINE exchange_x
      IMPLICIT NONE
      real(8), dimension(4, 3, -2:ny + 3):: send_buf_y, recv_buf_y

      IF (mod(position_x, 2) == 0) THEN
         mpitag = 1
         send_buf_y(:, 1, :) = u(:, 1, :)
         send_buf_y(:, 2, :) = u(:, 2, :)
         send_buf_y(:, 3, :) = u(:, 3, :)
         call MPI_send(send_buf_y, 12*(ny + 6), MPI_REAL8,&
        & destination(1), mpitag, MPI_COMM_WORLD, mpierror)

         mpitag = 1
         call MPI_RECV(recv_buf_y, 12*(ny + 6), MPI_REAL8, &
        & source(2), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, nx + 1, :) = recv_buf_y(:, 1, :)
         u(:, nx + 2, :) = recv_buf_y(:, 2, :)
         u(:, nx + 3, :) = recv_buf_y(:, 3, :)
      ELSE
         mpitag = 1
         call MPI_RECV(recv_buf_y, 12*(ny + 6), MPI_REAL8, &
        & source(2), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, nx + 1, :) = recv_buf_y(:, 1, :)
         u(:, nx + 2, :) = recv_buf_y(:, 2, :)
         u(:, nx + 3, :) = recv_buf_y(:, 3, :)
         mpitag = 1
         send_buf_y(:, 1, :) = u(:, 1, :)
         send_buf_y(:, 2, :) = u(:, 2, :)
         send_buf_y(:, 3, :) = u(:, 3, :)
         call MPI_send(send_buf_y, 12*(ny + 6), MPI_REAL8,&
        & destination(1), mpitag, MPI_COMM_WORLD, mpierror)
      END IF

      IF (mod(position_x, 2) == 0) THEN
         mpitag = 2
         send_buf_y(:, 1, :) = u(:, nx, :)
         send_buf_y(:, 2, :) = u(:, nx - 1, :)
         send_buf_y(:, 3, :) = u(:, nx - 2, :)
         call MPI_send(send_buf_y, 12*(ny + 6), MPI_REAL8,&
        & destination(2), mpitag, MPI_COMM_WORLD, mpierror)
         mpitag = 2
         call MPI_RECV(recv_buf_y, 12*(ny + 6), MPI_REAL8, &
        & source(1), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, 0, :) = recv_buf_y(:, 1, :)
         u(:, -1, :) = recv_buf_y(:, 2, :)
         u(:, -2, :) = recv_buf_y(:, 3, :)
      ELSE
         mpitag = 2
         call MPI_RECV(recv_buf_y, 12*(ny + 6), MPI_REAL8, &
        & source(1), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, 0, :) = recv_buf_y(:, 1, :)
         u(:, -1, :) = recv_buf_y(:, 2, :)
         u(:, -2, :) = recv_buf_y(:, 3, :)
         mpitag = 2
         send_buf_y(:, 1, :) = u(:, nx, :)
         send_buf_y(:, 2, :) = u(:, nx - 1, :)
         send_buf_y(:, 3, :) = u(:, nx - 2, :)
         call MPI_send(send_buf_y, 12*(ny + 6), MPI_REAL8,&
        & destination(2), mpitag, MPI_COMM_WORLD, mpierror)
      END IF
   END SUBROUTINE exchange_x

   SUBROUTINE exchange_y
      IMPLICIT NONE
      real(8), dimension(4, -2:nx + 3, 3):: send_buf_x, recv_buf_x

      IF (mod(position_y, 2) == 0) THEN
         mpitag = 4
         send_buf_x(:, :, 1) = u(:, :, 1)
         send_buf_x(:, :, 2) = u(:, :, 2)
         send_buf_x(:, :, 3) = u(:, :, 3)
         call MPI_send(send_buf_x, 12*(nx + 6), MPI_REAL8,&
        & destination(4), mpitag, MPI_COMM_WORLD, mpierror)
         mpitag = 4
         call MPI_RECV(recv_buf_x, 12*(nx + 6), MPI_REAL8, &
        & source(3), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, :, ny + 1) = recv_buf_x(:, :, 1)
         u(:, :, ny + 2) = recv_buf_x(:, :, 2)
         u(:, :, ny + 3) = recv_buf_x(:, :, 3)
      ELSE
         mpitag = 4
         call MPI_RECV(recv_buf_x, 12*(nx + 6), MPI_REAL8, &
        & source(3), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, :, ny + 1) = recv_buf_x(:, :, 1)
         u(:, :, ny + 2) = recv_buf_x(:, :, 2)
         u(:, :, ny + 3) = recv_buf_x(:, :, 3)
         mpitag = 4
         send_buf_x(:, :, 1) = u(:, :, 1)
         send_buf_x(:, :, 2) = u(:, :, 2)
         send_buf_x(:, :, 3) = u(:, :, 3)
         call MPI_send(send_buf_x, 12*(nx + 6), MPI_REAL8,&
        & destination(4), mpitag, MPI_COMM_WORLD, mpierror)
      END IF

      IF (mod(position_y, 2) == 0) THEN
         mpitag = 3
         send_buf_x(:, :, 1) = u(:, :, ny)
         send_buf_x(:, :, 2) = u(:, :, ny - 1)
         send_buf_x(:, :, 3) = u(:, :, ny - 2)
         call MPI_send(send_buf_x, 12*(nx + 6), MPI_REAL8,&
        & destination(3), mpitag, MPI_COMM_WORLD, mpierror)
         mpitag = 3
         call MPI_RECV(recv_buf_x, 12*(nx + 6), MPI_REAL8, &
        & source(4), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, :, 0) = recv_buf_x(:, :, 1)
         u(:, :, -1) = recv_buf_x(:, :, 2)
         u(:, :, -2) = recv_buf_x(:, :, 3)
      ELSE
         mpitag = 3
         call MPI_RECV(recv_buf_x, 12*(nx + 6), MPI_REAL8, &
        & source(4), mpitag, MPI_COMM_WORLD, send_status, mpierror)
         u(:, :, 0) = recv_buf_x(:, :, 1)
         u(:, :, -1) = recv_buf_x(:, :, 2)
         u(:, :, -2) = recv_buf_x(:, :, 3)
         mpitag = 3
         send_buf_x(:, :, 1) = u(:, :, ny)
         send_buf_x(:, :, 2) = u(:, :, ny - 1)
         send_buf_x(:, :, 3) = u(:, :, ny - 2)
         call MPI_send(send_buf_x, 12*(nx + 6), MPI_REAL8,&
        & destination(3), mpitag, MPI_COMM_WORLD, mpierror)
      END IF
   END SUBROUTINE exchange_y
   !**************************
   subroutine initdata
      implicit none
      real(8):: xleft, xright, yleft, yright, ue(4)
      real(8):: rho_L, v_L, p_L, rho_R, v_R, p_R
      REAL(8) :: x_center = 0.0, y_center = 0.0
      real(8) :: global_x, global_y
      real(8) :: total_dx = xm/total_x, total_dy = ym/total_y

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

      !Quarter domain:
      xleft = xm/cut_x*position_x
      xright = xm/cut_x*(position_x + 1)
      xlen = xright - xleft
      dx = xlen/nx
      yleft = ym/cut_y*position_y
      yright = ym/cut_y*(position_y + 1)
      ylen = yright - yleft
      dy = ylen/ny

      do i = -2, nx + 3
         x(i) = xleft + i*dx - dx/2.
      end do
      do i = -2, ny + 3
         y(i) = yleft + i*dy - dy/2.
      end do

      u = 0.
      u(1, :, :) = 1
      u(4, :, :) = 1d-12

      do j = 1, ny
      do i = 1, nx
          global_x = xleft + dx*(i - 0.5) !xL = 0
          global_y = yleft + dy*(j - 0.5) !yL = 0

          ! Check if this cell contains the blast center
          if (abs(global_x - x_center) <= dx/2.0 .AND. abs(global_y - y_center) <= dy/2.0) then
              ! Found the blast center cell, inject energy here
              u(4, i, j) = 0.244816 / (dx * dy)

          end if
      end do
  end do

      call bc
   end subroutine initdata
   !******************************

   !******************************
   SUBROUTINE output
      IMPLICIT NONE
      ! The global grid dimensions should be defined here or passed as arguments
      !INTEGER, PARAMETER :: total_x = 800, total_y = 400
      
      ! Define  using the new global dimensions
      real(8):: er1, eri, er1old, eriold, o1, oi, er1temp, eritemp
     ! REAL(8):: utemp(4, 1:total_x, 1:total_y)
      
      REAL(8):: recv_buf(4, 1:nx, 1:ny)
      REAL(8), ALLOCATABLE :: buffer(:)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: buffer_addr
      INTEGER :: buffer_size_bytes
    !  INTEGER :: mpierror

     real(8):: lrho, lm, ln, le, pres
      integer:: den_tmp, pre_tmp

! Calculate the size of the buffer in bytes
    ! (size of REAL(8) is 8 bytes)
    buffer_size_bytes = (4 * (cut_x * nx + 6) * (cut_y * ny + 6)) * SIZEOF(1.0_8)

! Allocate the buffer
    ALLOCATE(buffer(buffer_size_bytes))

      er1 = 0.
      eri = 0.
  IF (KIND(den_count) /= KIND(den_tmp)) THEN
     PRINT *, 'Error: den_count and den_tmp kinds do not match!'
     PRINT *, '  den_count kind =', KIND(den_count)
     PRINT *, '  den_tmp kind =', KIND(den_tmp)
     STOP
  END IF
      call MPI_AllReduce(den_count, den_tmp, 1, MPI_INTEGER, MPI_SUM, &
      & MPI_COMM_WORLD, mpierror)

      call MPI_AllReduce(pre_count, pre_tmp, 1, MPI_INTEGER, MPI_SUM, &
      & MPI_COMM_WORLD, mpierror)
      den_count = den_tmp
      pre_count = pre_tmp

     ! call MPI_BUFFER_ATTACH(buffer, 4*total*total, mpierror)
   ! CALL MPI_BUFFER_ATTACH(buffer, buffer_size_bytes, mpierror)
     call gather_to_root

 !if (mpirank .ne. 0) then
    !     mpitag = mpirank
    !     mpides = 0
    !     call MPI_send(u(:, 1:nx, 1:ny), 4*nx*ny, MPI_REAL8, mpides,&
    !            & mpitag, MPI_COMM_WORLD, mpierror)
   !   else
   !      do ig = 1, nx
   !      do jg = 1, ny
   !         utemp(:, ig, jg) = u(:, ig, jg)
   !      end do
   !      end do
   !      DO i = 1, mpisize - 1
   !         mpisource = i
   !         mpitag = i
   !         call MPI_RECV(recv_buf, 4*nx*ny, &
   !                        MPI_REAL8, mpisource, mpitag, MPI_COMM_WORLD, send_status, mpierror)
   !         DO ig = 1, nx
   !         DO jg = 1, ny
   !            utemp(:, ig + mod(i, cut_x)*nx, jg + i/cut_x*ny) = recv_buf(:, ig, jg)
   !         END DO
   !         END DO
   !      END DO
   !   end if
!      call MPI_BUFFER_DETACH(buffer, 4*total*total, mpierror)
   ! CALL MPI_BUFFER_DETACH(buffer, buffer_size_bytes, mpierror)
   ! DEALLOCATE(buffer)

      if (mpirank == 0) then

         ! open(73, file='den.dat')
         ! open(74, file='pre.dat')
         ! open(75, file='temp.dat')
         !Create new file paths for stored variables
         write(den_file_path, '(A,A)') trim(output_dir), "/den.dat"
         write(m1_file_path, '(A,A)') trim(output_dir), "/m1.dat"
         write(m2_file_path, '(A,A)') trim(output_dir), "/m2.dat"
         write(pre_file_path, '(A,A)') trim(output_dir), "/pre.dat"

         open (73, file=trim(den_file_path))
         open(74, file=trim(pre_file_path))
         open(75, file=trim(m1_file_path))
         open(76, file=trim(m2_file_path))   


         write (73, *) 'VARIABLES = "X", "Y", "DENSITY"'

         write (73, *) 'zone i=', total_x, ', j=', total_y
         write (74, *) 'zone i=', total_x, ', j=', total_y
         write (75, *) 'zone i=', total_x, ', j=', total_y
         write (76, *) 'zone i=', total_x, ', j=', total_y

        ! do j = 1, total_y
        ! do i = 1, total_x

        !    lrho = utemp(1, i, total_y - j + 1)

        !    lm = utemp(2, i, total_y - j + 1)
        !    ln = utemp(3, i, total_y - j + 1)
        !    le = utemp(4, i, total_y - j + 1)
        !    pres = gamma1*(le - 0.5*lm**2/lrho - 0.5*ln**2/lrho)

        !    write (73, *) dx*(i - 0.5), dy*(j - 0.5) - 0.25, log(abs(lrho))
        !    write (74, *) dx*(i - 0.5), dy*(j - 0.5) - 0.25, log(abs(pres))
        !    write (75, *) dx*(i - 0.5), dy*(j - 0.5) - 0.25, log(abs(pres/lrho))
        ! end do
        ! end do

         do j = 1, total_y
         do i = 1, total_x

            lrho = utemp(1, i, j)

            lm = utemp(2, i, j)
            ln = utemp(3, i, j)
            le = utemp(4, i, j)
            pres = gamma1*(le - 0.5*lm**2/lrho - 0.5*ln**2/lrho)

            write (73, *) dx*(i - 0.5), dy*(j - 0.5), lrho !Shift to cell center (finite volume)
            write (74, *) dx*(i - 0.5), dy*(j - 0.5), pres
            write (75, *) dx*(i - 0.5), dy*(j - 0.5), lm
            write (76, *) dx*(i - 0.5), dy*(j - 0.5), ln
         end do
         end do
         close (73)
         close (74)
         close (75)
         close (76)

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
   FUNCTION nod_lim(v_1,v_2,v_bar)
      IMPLICIT NONE
      REAL(8), DIMENSION(4), INTENT(IN):: v_1,v_2,v_bar
      REAL(8):: p_v_1, p_v_2, p_v_bar, t,t_1, delta, d_1_2, d_1_bar
      REAL(8), DIMENSION(4):: nod_lim


    P_v_1 = (gamma1)*(v_1(4) - 0.5*(v_1(2)*v_1(2) + v_1(3)*v_1(3))/v_1(1));
    P_v_2 = (gamma1)*(v_2(4) - 0.5*(v_2(2)*v_2(2) + v_2(3)*v_2(3))/v_2(1));
   

    !if (p_v_2 > 0 .AND. v_2(1)>0) then
    if (p_v_2 > eps) then
        t = (P_v_1 - eps)/(P_v_1 - P_v_2);
        t_pos_count = t_pos_count + 1 

    else
        P_v_bar = (gamma1)*(v_bar(4) - 0.5*(v_bar(2)*v_bar(2) + v_bar(3)*v_bar(3))/v_bar(1));
       ! t_1 = P_v_1/(P_v_1 - P_v_bar+1e-15);
        t_1 = (P_v_1-eps)/(P_v_1 - P_v_bar);
        d_1_2 = sqrt((v_2(4) - v_1(4))**2 + (v_2(3) - v_1(3))**2 + (v_2(2) - v_1(2))**2 + (v_2(1) - v_1(1))**2);
        d_1_bar = sqrt((v_bar(4) - v_1(4))**2 + (v_bar(3) - v_1(3))**2 + (v_bar(2) - v_1(2))**2 + (v_bar(1) - v_1(1))**2);
        
        t = min(t_1*d_1_bar/d_1_2,0.25);

         t_neg_count = t_neg_count + 1

    
    end if

    nod_lim(1) = (1 - t)*v_1(1) + t*v_2(1);
    nod_lim(2) = (1 - t)*v_1(2) + t*v_2(2);
    nod_lim(3) = (1 - t)*v_1(3) + t*v_2(3);
    nod_lim(4) = (1 - t)*v_1(4) + t*v_2(4);


      RETURN
   END FUNCTION nod_lim


   !*****************************
    !*****************************
  SUBROUTINE sweep_lim
  IMPLICIT NONE
  REAL(8), dimension(4)::  v_j, v_i, v_p, u_bar, v_j_t, v_i_t, v_p_t, d
  REAL(8), dimension(1:total_x, 1:total_y):: rho_t, m1_t, m2_t, E_t, P_check
  REAL(8):: rho_avg, m1_avg, m2_avg, E_avg, pres_v_j , w_p, w_j, pres_v_i , w_i, rho_min, p_min
  INTEGER:: q

  REAL(8)::  pres_v_p,  pres_v_i_t, pres_v_p_t, dist_tilde, dist_pi
 
if (mpirank == 0) then  
 
  call scalar_sweep


  DO j = 1, total_y
     DO i = 1, total_x
           rho_t(i,j) = utemp(1,i,j);
           m1_t(i,j)  = utemp(2,i,j);
           m2_t(i,j)  = utemp(3,i,j);
           E_t(i,j)   = utemp(4,i,j);
     END DO
  END DO


  rho_avg = 0.0
  m1_avg = 0.0
  m2_avg = 0.0
  E_avg = 0.0


  do j = 1, total_y
     do i = 1, total_x
        rho_avg = rho_avg + weights(i,j) * rho_t(i,j)
        m1_avg = m1_avg + weights(i,j) * m1_t(i,j)
        m2_avg = m2_avg + weights(i,j) * m2_t(i,j)
        E_avg = E_avg + weights(i,j) * E_t(i,j)
     end do
  end do


  u_bar(1) = rho_avg
  u_bar(2) = m1_avg
  u_bar(3) = m2_avg
  u_bar(4) = E_avg


if (gamma1*(u_bar(4) - 0.5*(u_bar(2)*u_bar(2) + u_bar(3)*u_bar(3))/u_bar(1)) <= 0.) then
    if (mpirank == 0) then 
    print *, "ERROR: P(u_bar) <= 0."
    print *, "P(u_bar) = ", (gamma1*(u_bar(4) - 0.5*(u_bar(2)*u_bar(2) + u_bar(3)*u_bar(3))/u_bar(1)))
    end if 

    stop "Program terminated due to incorrect physical conditions."
end if

   sweep_count = sweep_count + 1
   t_pos_count = 0
   t_neg_count = 0

   DO q = 1, swp_num
     if (mod(q,2) == 1) then  !Alternate sweep 1 and sweep 2
        !call scalar_sweep


      
      
           !Forward Sweep
           DO i = 1, (total_x*total_y-1)
                 v_i(1) = utemp(1,sweepI_i(i),sweepI_j(i));
                 v_i(2) = utemp(2,sweepI_i(i),sweepI_j(i));
                 v_i(3) = utemp(3,sweepI_i(i),sweepI_j(i));
                 v_i(4) = utemp(4,sweepI_i(i),sweepI_j(i));
    
                 pres_v_i = gamma1*(v_i(4) - 0.5*(v_i(2)*v_i(2) + v_i(3)*v_i(3))/(v_i(1)))
                 w_i = weights(sweepI_i(i),sweepI_j(i));
    
                 !if (pres_v_i <= 0. .OR. v_i(1) <= 0.) then
                 if (pres_v_i <= eps .OR. v_i(1) <= eps) then

                       v_p(1) = utemp(1,sweepI_i(i+1),sweepI_j(i+1));
                       v_p(2) = utemp(2,sweepI_i(i+1),sweepI_j(i+1));
                       v_p(3) = utemp(3,sweepI_i(i+1),sweepI_j(i+1));
                       v_p(4) = utemp(4,sweepI_i(i+1),sweepI_j(i+1));
                      
                       w_p = weights(sweepI_i(i+1),sweepI_j(i+1));
                       v_i_t = nod_lim(v_i,v_p,u_bar);
    
                          d  = v_i - v_i_t;
                          v_p_t =  v_p  + (w_i/w_p)*(d);
                         
                          utemp(1,sweepI_i(i),sweepI_j(i)) = v_i_t(1);
                          utemp(2,sweepI_i(i),sweepI_j(i)) = v_i_t(2);
                          utemp(3,sweepI_i(i),sweepI_j(i)) = v_i_t(3);
                          utemp(4,sweepI_i(i),sweepI_j(i)) = v_i_t(4);
    
                          utemp(1,sweepI_i(i+1),sweepI_j(i+1)) = v_p_t(1);
                          utemp(2,sweepI_i(i+1),sweepI_j(i+1)) = v_p_t(2);
                          utemp(3,sweepI_i(i+1),sweepI_j(i+1)) = v_p_t(3);
                          utemp(4,sweepI_i(i+1),sweepI_j(i+1)) = v_p_t(4);


                          pres_v_p = gamma1*(v_p(4) - 0.5*(v_p(2)*v_p(2) + v_p(3)*v_p(3))/v_p(1))
                          pres_v_i_t = gamma1*(v_i_t(4) - 0.5*(v_i_t(2)*v_i_t(2) + v_i_t(3)*v_i_t(3))/v_i_t(1))
                          pres_v_p_t = gamma1*(v_p_t(4) - 0.5*(v_p_t(2)*v_p_t(2) + v_p_t(3)*v_p_t(3))/v_p_t(1));



                         !   if (q == 1) then
                         !  end if

                          ! if (q == 1) then
                          ! end if

                        ! if (pres_v_i + pres_v_p > pres_v_i_t + pres_v_p_t + eps) then

                           ! Distance between v_p and v_i
                          !dist_pi = sqrt(sum((v_p - v_i)**2))

                            ! Distance between v_p_t and v_i_t
                          !dist_tilde = sqrt(sum((v_p_t - v_i_t)**2))

                             ! stop
                         !end if 
                          ! if (i == 66 .AND. y_iter == 64) then
                          !end if
                 end if
           end do
       
           !Backward Sweep
           DO i = total_x*total_y, 2, -1
          
                 v_i(1) = utemp(1,sweepI_i(i),sweepI_j(i));
                 v_i(2) = utemp(2,sweepI_i(i),sweepI_j(i));
                 v_i(3) = utemp(3,sweepI_i(i),sweepI_j(i));
                 v_i(4) = utemp(4,sweepI_i(i),sweepI_j(i));
    
                 pres_v_i =  gamma1*(v_i(4) - 0.5*(v_i(2)*v_i(2) + v_i(3)*v_i(3))/v_i(1));




                 w_i = weights(sweepI_i(i),sweepI_j(i));
    
                 !if (pres_v_i <= 0.  .OR. v_i(1) <= 0.) then
                 if (pres_v_i <= eps  .OR. v_i(1) <= eps) then

                    v_p(1) = utemp(1,sweepI_i(i-1),sweepI_j(i-1));
                    v_p(2) = utemp(2,sweepI_i(i-1),sweepI_j(i-1));
                    v_p(3) = utemp(3,sweepI_i(i-1),sweepI_j(i-1));
                    v_p(4) = utemp(4,sweepI_i(i-1),sweepI_j(i-1));
                   
                    w_p = weights(sweepI_i(i-1),sweepI_j(i-1));
                    v_i_t = nod_lim(v_i,v_p,u_bar);
    
                          d  = v_i - v_i_t;
    
                          v_p_t =  v_p  + (w_i/ w_p)*(d);
                         
                          utemp(1,sweepI_i(i),sweepI_j(i)) = v_i_t(1);
                          utemp(2,sweepI_i(i),sweepI_j(i)) = v_i_t(2);
                          utemp(3,sweepI_i(i),sweepI_j(i)) = v_i_t(3);
                          utemp(4,sweepI_i(i),sweepI_j(i)) = v_i_t(4);
    
                          utemp(1,sweepI_i(i-1),sweepI_j(i-1)) = v_p_t(1);
                          utemp(2,sweepI_i(i-1),sweepI_j(i-1)) = v_p_t(2);
                          utemp(3,sweepI_i(i-1),sweepI_j(i-1)) = v_p_t(3);
                          utemp(4,sweepI_i(i-1),sweepI_j(i-1)) = v_p_t(4);


                          pres_v_p = gamma1*(v_p(4) - 0.5*(v_p(2)*v_p(2) + v_p(3)*v_p(3))/v_p(1))
                          pres_v_i_t = gamma1*(v_i_t(4) - 0.5*(v_i_t(2)*v_i_t(2) + v_i_t(3)*v_i_t(3))/v_i_t(1))
                          pres_v_p_t = gamma1*(v_p_t(4) - 0.5*(v_p_t(2)*v_p_t(2) + v_p_t(3)*v_p_t(3))/v_p_t(1));

                       !  if (pres_v_i + pres_v_p > pres_v_i_t + pres_v_p_t + eps) then
                                 ! Distance between v_p and v_i
                        !      dist_pi = sqrt(sum((v_p - v_i)**2))

                                 ! Distance between v_p_t and v_i_t
                        !      dist_tilde = sqrt(sum((v_p_t - v_i_t)**2))

                           !   stop
                             ! stop
                        ! end if 
    
                 end if
           end do




     else


    
         !Forward Sweep
          DO i = 1, (total_x*total_y-1)
  
                v_i(1) = utemp(1,sweepII_i(i),sweepII_j(i));
                v_i(2) = utemp(2,sweepII_i(i),sweepII_j(i));
                v_i(3) = utemp(3,sweepII_i(i),sweepII_j(i));
                v_i(4) = utemp(4,sweepII_i(i),sweepII_j(i));
  
                pres_v_i = gamma1*(v_i(4) - 0.5*(v_i(2)*v_i(2) + v_i(3)*v_i(3))/(v_i(1)))
                w_i = weights(sweepII_i(i),sweepII_j(i));
  
                !if (pres_v_i <= 0.0 .OR. v_i(1) <= 0.0) then
                if (pres_v_i <= eps .OR. v_i(1) <= eps) then

                      v_p(1) = utemp(1,sweepII_i(i+1),sweepII_j(i+1));
                      v_p(2) = utemp(2,sweepII_i(i+1),sweepII_j(i+1));
                      v_p(3) = utemp(3,sweepII_i(i+1),sweepII_j(i+1));
                      v_p(4) = utemp(4,sweepII_i(i+1),sweepII_j(i+1));
                    
                      w_p = weights(sweepII_i(i+1),sweepII_j(i+1));
                      v_i_t = nod_lim(v_i,v_p,u_bar);
  
                         d  = v_i - v_i_t;
                         v_p_t =  v_p  + (w_i/w_p)*(d);
                       
                         utemp(1,sweepII_i(i),sweepII_j(i)) = v_i_t(1);
                         utemp(2,sweepII_i(i),sweepII_j(i)) = v_i_t(2);
                         utemp(3,sweepII_i(i),sweepII_j(i)) = v_i_t(3);
                         utemp(4,sweepII_i(i),sweepII_j(i)) = v_i_t(4);
  
                         utemp(1,sweepII_i(i+1),sweepII_j(i+1)) = v_p_t(1);
                         utemp(2,sweepII_i(i+1),sweepII_j(i+1)) = v_p_t(2);
                         utemp(3,sweepII_i(i+1),sweepII_j(i+1)) = v_p_t(3);
                         utemp(4,sweepII_i(i+1),sweepII_j(i+1)) = v_p_t(4);




                          pres_v_p = gamma1*(v_p(4) - 0.5*(v_p(2)*v_p(2) + v_p(3)*v_p(3))/v_p(1))
                          pres_v_i_t = gamma1*(v_i_t(4) - 0.5*(v_i_t(2)*v_i_t(2) + v_i_t(3)*v_i_t(3))/v_i_t(1))
                          pres_v_p_t = gamma1*(v_p_t(4) - 0.5*(v_p_t(2)*v_p_t(2) + v_p_t(3)*v_p_t(3))/v_p_t(1));
                          
                    !      if (pres_v_i + pres_v_p > pres_v_i_t + pres_v_p_t+eps) then


                           ! Distance between v_p and v_i
                    !      dist_pi = sqrt(sum((v_p - v_i)**2))

                            ! Distance between v_p_t and v_i_t
                     !     dist_tilde = sqrt(sum((v_p_t - v_i_t)**2))

                            !  stop
                     !    end if 




                         ! if (i == 66 .AND. y_iter == 64) then
                         !end if
                end if
          end do
     
          !Backward Sweep
          DO i = (total_x*total_y), 2, -1
        
                v_i(1) = utemp(1,sweepII_i(i),sweepII_j(i))
                v_i(2) = utemp(2,sweepII_i(i),sweepII_j(i))
                v_i(3) = utemp(3,sweepII_i(i),sweepII_j(i))
                v_i(4) = utemp(4,sweepII_i(i),sweepII_j(i))
  
                pres_v_i =  gamma1*(v_i(4) - 0.5*(v_i(2)*v_i(2) + v_i(3)*v_i(3))/v_i(1))

                w_i = weights(sweepII_i(i),sweepII_j(i));
  
                if (pres_v_i <= eps  .OR. v_i(1) <= eps) then
  
                   v_p(1) = utemp(1,sweepII_i(i-1),sweepII_j(i-1));
                   v_p(2) = utemp(2,sweepII_i(i-1),sweepII_j(i-1));
                   v_p(3) = utemp(3,sweepII_i(i-1),sweepII_j(i-1));
                   v_p(4) = utemp(4,sweepII_i(i-1),sweepII_j(i-1));
                 
                   w_p = weights(sweepII_i(i-1),sweepII_j(i-1));
                   v_i_t = nod_lim(v_i,v_p,u_bar);
  
                         d  = v_i - v_i_t;
  
                         v_p_t =  v_p  + (w_i/ w_p)*(d);
                       
                         utemp(1,sweepII_i(i),sweepII_j(i)) = v_i_t(1);
                         utemp(2,sweepII_i(i),sweepII_j(i)) = v_i_t(2);
                         utemp(3,sweepII_i(i),sweepII_j(i)) = v_i_t(3);
                         utemp(4,sweepII_i(i),sweepII_j(i)) = v_i_t(4);
  
                         utemp(1,sweepII_i(i-1),sweepII_j(i-1)) = v_p_t(1);
                         utemp(2,sweepII_i(i-1),sweepII_j(i-1)) = v_p_t(2);
                         utemp(3,sweepII_i(i-1),sweepII_j(i-1)) = v_p_t(3);
                         utemp(4,sweepII_i(i-1),sweepII_j(i-1)) = v_p_t(4);


                          pres_v_p = gamma1*(v_p(4) - 0.5*(v_p(2)*v_p(2) + v_p(3)*v_p(3))/v_p(1))
                          pres_v_i_t = gamma1*(v_i_t(4) - 0.5*(v_i_t(2)*v_i_t(2) + v_i_t(3)*v_i_t(3))/v_i_t(1))
                          pres_v_p_t = gamma1*(v_p_t(4) - 0.5*(v_p_t(2)*v_p_t(2) + v_p_t(3)*v_p_t(3))/v_p_t(1));

                     !    if (pres_v_i + pres_v_p > pres_v_i_t + pres_v_p_t + eps) then


                           ! Distance between v_p and v_i
                      !    dist_pi = sqrt(sum((v_p - v_i)**2))

                            ! Distance between v_p_t and v_i_t
                      !    dist_tilde = sqrt(sum((v_p_t - v_i_t)**2))

                            !  stop
                      !   end if 
  
                end if
          end do
      end if

      rho_min = 1d20
      P_min = 1d20
      DO j = 1, total_y
         DO i = 1, total_x
            !P_check(j) = gamma1*(utemp(4,x_iter,j) - 0.5*(utemp(2,x_iter,j)**2 + utemp(3,x_iter,j)**2)/utemp(1,x_iter,j));
            P_min = min(gamma1*(utemp(4,i,j) - 0.5*(utemp(2,i,j)*utemp(2,i,j) + utemp(3,i,j)*utemp(3,i,j))/utemp(1,i,j)),P_min);
            rho_min = min(utemp(1,i,j), rho_min);
         
         end do
      end do

      pre_count = pre_count + 1 
      


     
     ! if (mpirank == 0) then 
     ! end if 

     !if (q == 10) then
     ! stop
     ! end if

      if (P_min > eps .AND. rho_min>eps) then !(P_min > 0 .AND. rho_min>0) then

          sweep_track(sweep_count) = q
            neg_track(sweep_count) = t_neg_count
            pos_track(sweep_count) = t_pos_count

         exit
      end if


      if (q == swp_num) then
         if (mpirank == 0) then
         print *, "SWEEP FAIL: swp_num =",swp_num
         end if
         stop
      end if


     end do
  end if
     
  postag = 0


     END SUBROUTINE  sweep_lim


!*****************************
      SUBROUTINE scalar_sweep
         IMPLICIT NONE
         REAL(8)::  v_j, v_jp1, v_jm1,  v_i, v_ip1, v_im1
         REAL(8), dimension(1:total_x,1:total_y):: rho_t 
         REAL(8)::  w_j, w_jp1, w_jm1, w_i, w_ip1, w_im1, delta, rho_min, rho_max
         INTEGER:: ss
               
      if (use_limiter == 2) then
      rho_min = 1d20
      DO j = 1, total_y
         DO i = 1, total_x
             rho_min = min(utemp(1,i,j), rho_min);
         end do
      end do

      delta = 1d-12;
      if (rho_min < delta) then
      DO j = 1, total_y
         DO i = 1, total_x
            rho_t(i,j) = utemp(1,i,j);
         end do
      end do
      
      DO i = 1, (total_x*total_y-1)
         v_i = rho_t(sweepI_i(i),sweepI_j(i)); v_ip1 = rho_t(sweepI_i(i+1),sweepI_j(i+1)); w_i = weights(sweepI_i(i),sweepI_j(i)); w_ip1 = weights(sweepI_i(i+1),sweepI_j(i+1));
         if (v_i < delta) then
            rho_t(sweepI_i(i+1),sweepI_j(i+1)) =  v_ip1  + (w_i/w_ip1)*(v_i - delta);
            rho_t(sweepI_i(i),sweepI_j(i)) = delta;
            end if
      end do


      DO i = total_x*total_y, 2, -1
         v_i = rho_t(sweepI_i(i),sweepI_j(i)); v_im1 = rho_t(sweepI_i(i-1),sweepI_j(i-1)); w_i = weights(sweepI_i(i),sweepI_j(i)); w_im1 = weights(sweepI_i(i-1),sweepI_j(i-1));
         if (v_i < delta) then
            rho_t(sweepI_i(i-1),sweepI_j(i-1)) =  v_im1  + (w_i/ w_im1)*(v_i - delta);
            rho_t(sweepI_i(i),sweepI_j(i)) = delta;
         end if
      end do

      !Fill in values for u!
      DO  j = 1, total_y
         DO i = 1, total_x
            utemp(1,i,j) = rho_t(i,j);
         end do
      end do

      den_count = den_count + 1;
   end if 
      rho_min = 1d20
      DO j = 1, total_y
         DO i = 1, total_x
            rho_min = min(utemp(1,i,j), rho_min);
         end do
      end do

      rho_max = 0
      DO j = 1, total_y
         DO i = 1, total_x
            rho_max = max(utemp(1,i,j), rho_max);
         end do
      end do

      postag = 0

   end if
            
      
 END SUBROUTINE  scalar_sweep  


SUBROUTINE sweeping_order
IMPLICIT NONE


! ====================
! Sweep I: column-by-column
! ====================
k = 0
do j = 1, total_y
if (mod(j,2) == 1) then
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
if (mod(i,2) == 1) then
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
   REAL(8):: recv_buf(4, 1:nx, 1:ny)
   REAL(8), ALLOCATABLE :: buffer(:)
   INTEGER(KIND=MPI_ADDRESS_KIND) :: buffer_addr
   INTEGER :: buffer_size_bytes
   !INTEGER :: mpierror

   real(8):: lrho, lm, ln, le, pres
   integer:: den_tmp, pre_tmp

! Calculate the size of the buffer in bytes
    ! (size of REAL(8) is 8 bytes)
    buffer_size_bytes = (4 * (cut_x * nx + 6) * (cut_y * ny + 6)) * SIZEOF(1.0_8)

! Allocate the buffer
    ALLOCATE(buffer(buffer_size_bytes))

   ! Assumes these variables are already declared in the program:
   ! u, utemp, nx, ny, cut_x, cut_y, mpirank, mpisize, recv_buf, send_status
   ! buffer, buffer_size_bytes, mpierror, mpitag, mpides, mpisource

   ! Attach buffer for buffered sends
   CALL MPI_BUFFER_ATTACH(buffer, buffer_size_bytes, mpierror)

   if (mpirank .ne. 0) then
      mpitag = mpirank
     mpides = 0
      CALL MPI_SEND(u(:,1:nx,1:ny), 4*nx*ny, MPI_REAL8, mpides, mpitag, MPI_COMM_WORLD, mpierror)
   else
   !   ! Root copies its own block
      DO ig = 1, nx
         DO jg = 1, ny
            utemp(:, ig, jg) = u(:, ig, jg)
         END DO
      END DO

      ! Receive from all other ranks
      DO i = 1, mpisize - 1
         mpisource = i
         mpitag    = i
         CALL MPI_RECV(recv_buf, 4*nx*ny, MPI_REAL8, mpisource, mpitag, MPI_COMM_WORLD, send_status, mpierror)!
         DO ig = 1, nx
            DO jg = 1, ny
               utemp(:, ig + mod(i, cut_x)*nx, jg + i/cut_x*ny) = recv_buf(:, ig, jg)
            END DO
         END DO
      END DO
   END IF

 !  if (mpirank == 0) then
!end if

 !  if (mpirank == 0) then
!    ! Copy own block
!    do ig = 1, nx
!        do jg = 1, ny
!            utemp(:, ig, jg) = u(:, ig, jg)
!        end do
!    end do

    ! Receive blocks from other ranks
!    do i = 1, mpisize-1
!        mpisource = i
 !       mpitag = i
 !       call MPI_RECV(recv_buf, 4*nx*ny, MPI_REAL8, mpisource, mpitag, MPI_COMM_WORLD, send_status, mpierror)

!        do ig = 1, nx
!            do jg = 1, ny
!                utemp(:, ig + mod(i, cut_x)*nx, jg + i/cut_x*ny) = recv_buf(:, ig, jg)
!            end do
!        end do
!    end do
!end if

   ! Detach and deallocate buffer
   CALL MPI_BUFFER_DETACH(buffer, buffer_size_bytes, mpierror)
   DEALLOCATE(buffer)
end subroutine gather_to_root


!-------------------------------------------------------------
subroutine scatter_from_root
   implicit none
   REAL(8):: recv_buf(4, 1:nx, 1:ny)
   REAL(8) :: send_buf(4, 1:nx, 1:ny)

      REAL(8), ALLOCATABLE :: buffer(:)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: buffer_addr
      INTEGER :: buffer_size_bytes
      !INTEGER :: mpierror

     real(8):: lrho, lm, ln, le, pres
      integer:: den_tmp, pre_tmp

! Calculate the size of the buffer in bytes
    ! (size of REAL(8) is 8 bytes)
    buffer_size_bytes = (4 * (cut_x * nx + 6) * (cut_y * ny + 6)) * SIZEOF(1.0_8)

! Allocate the buffer
    ALLOCATE(buffer(buffer_size_bytes))

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
         mpitag  = i

         DO ig = 1, nx
            DO jg = 1, ny
               send_buf(:, ig, jg) = utemp(:, ig + mod(i, cut_x)*nx, jg + i/cut_x*ny)
            END DO
         END DO

         CALL MPI_BSEND(send_buf, 4*nx*ny, MPI_REAL8, mpides, mpitag, MPI_COMM_WORLD, mpierror)
      END DO

   else
      ! Non-root receives its block
      mpisource = 0
      mpitag    = mpirank
      CALL MPI_RECV(recv_buf, 4*nx*ny, MPI_REAL8, mpisource, mpitag, MPI_COMM_WORLD, recv_status, mpierror)

      DO ig = 1, nx
         DO jg = 1, ny
            u(:, ig, jg) = recv_buf(:, ig, jg)
         END DO
      END DO
   END IF

   ! Detach and deallocate buffer
   CALL MPI_BUFFER_DETACH(buffer, buffer_size_bytes, mpierror)
   DEALLOCATE(buffer)
end subroutine scatter_from_root



END PROGRAM FDWENO4EULER1D