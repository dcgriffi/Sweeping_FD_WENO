PROGRAM FDWENO4EULER1D
   USE weno5
   IMPLICIT NONE !80, 160, 320, 640
   INTEGER, PARAMETER:: nx= 800 , ny = 200, n = nx, ierror = 1, use_limiter = 2, swp_num = 20000
!time_terminal = 0.001
   REAL(8), PARAMETER::  time_terminal = 0.001, gamma = 5./3., gamma1 = gamma - 1., eps = 1d-13,eps2=1d-8, xm=1., ym=0.25, ujet=800.
   REAL(8), DIMENSION(4, -2:nx + 3, -2:ny + 3):: u, u0, RHS, fp_flux, fn_flux, fp, fn, gp_flux, gn_flux, gp, gn
   REAL(8), DIMENSION(-2:nx + 3):: x
   REAL(8), DIMENSION(-2:ny + 3):: y
   integer:: rktag, kcount
   REAL(8), DIMENSION(-2:nx + 3, -2:ny + 3)::alpha
   real(8):: dx, dy, dt, time, alphamax
   integer:: ig, i, j, k, kk, m, den_count = 0, pre_count = 0, postag, x_iter, y_iter

   !Added for 2D sweeping
   INTEGER, PARAMETER:: total_x = nx, total_y = ny
   Real(8), PARAMETER:: nSteps_x = total_x, nsteps_y = total_y
   REAL(8), DIMENSION(1:total_x, 1:total_y):: weights
   integer, allocatable :: sweepI_i(:), sweepI_j(:)
   integer, allocatable :: sweepII_i(:), sweepII_j(:)

   !Adding sweep data
   INTEGER, PARAMETER :: total_max = 20000
   INTEGER, DIMENSION(1:total_max) :: sweep_track = 0
   INTEGER :: sweep_count = 0, total_sweep = 0, max_sweep = 0
   real(8):: avg_sweep

   INTEGER, DIMENSION(1:total_max) :: neg_track = 0
   INTEGER, DIMENSION(1:total_max) :: pos_track = 0
   INTEGER :: t_neg_count, t_pos_count, total_neg, total_pos, max_neg, max_pos
   real(8):: avg_pos, avg_neg



   ! Allocate arrays to store indices in sweep order (do in Rank 1 MPI)
   allocate(sweepI_i(total_x*total_y), sweepI_j(total_x*total_y))
   allocate(sweepII_i(total_x*total_y), sweepII_j(total_x*total_y))
   call sweeping_order

   call initdata

   DO i = 1, total_x 
      DO j = 1, total_y
      weights(i,j) = 1./(nSteps_x*nSteps_y)
      END DO
   END DO

   time = 0d0
   kcount = 0

   DO WHILE (time <= time_terminal - 1d-10)
      CALL rk
      time = time + dt
      kcount = kcount + 1
      if (kcount/20*20 .eq. kcount) write (6, *) kcount, time/time_terminal*100d0, '%'
   END DO
   CALL output

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

      if (rho < 1d-8) vx = 0.
      if (rho < 1d-8) vy = 0.

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

      if (rho < 1d-8) vx = 0.
      if (rho < 1d-8) vy = 0.

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
      REAL(8):: den, xmt, ymt, eng, t0, t1, t2, t3, pre, vex, vey, vxm, vym, hm, qm, cm, rcm, b1, b2
      REAL(8):: ul(4), ur(4), d_ave, m_ave, e_ave, p_ave, theta, point(4, 2), d_l, m_l, e_l, p_l, n_l, n_ave
      call findalpha

      do i = -2, nx + 3
      do j = -2, ny + 3

         fp(:, i, j) = 0.5*(u(:, i, j) + fflux(u(:, i, j))/alphamax)
         fn(:, i, j) = 0.5*(u(:, i, j) - fflux(u(:, i, j))/alphamax)
         gp(:, i, j) = 0.5*(u(:, i, j) + gflux(u(:, i, j))/alphamax)
         gn(:, i, j) = 0.5*(u(:, i, j) - gflux(u(:, i, j))/alphamax)
      end do
      end do

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
         DO i = 0, nx + 1
            !********************************
            ! Roe Average Matrix for f'(u)
            !********************************

            t0 = w(i)/(w(i) + w(i + 1) + 1.e-15)
            t1 = 1.-t0
            vxm = t0*vx(i) + t1*vx(i + 1)
            vym = t0*vy(i) + t1*vy(i + 1)
            hm = t0*h(i) + t1*h(i + 1)
            qm = 0.5*(vxm*vxm + vym*vym)
            cm = sqrt(abs(gamma1*(hm - qm))) + 1.e-15
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

!          evr=0.
!          evl=0.
!          do ig=1,4
!             evr(ig,ig)=1.
!             evl(ig,ig)=1.
!          enddo
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
            !*****************************
            ! positivity limiter
            !*****************************
            !if (fp(1,i,j)<0. .or. fp(4,i,j)-fp(2,i,j)**2/fp(1,i,j)/2.-fp(3,i,j)**2/fp(1,i,j)/2.<0.) then
            !else
            !endif

            !if (fn(1,i,j)<0. .or. fn(4,i,j)-fn(2,i,j)**2/fn(1,i,j)/2.-fn(3,i,j)**2/fn(1,i,j)/2.<0.) then
            !else
            !endif

			if (use_limiter == 1) then

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

            if (theta < 1d0) then
               pre_count = pre_count + 1
               fn_flux(1, i, j) = theta*fn_flux(1, i, j) + (1.-theta)*d_ave
               fn_flux(2, i, j) = theta*fn_flux(2, i, j) + (1.-theta)*m_ave
               fn_flux(3, i, j) = theta*fn_flux(3, i, j) + (1.-theta)*n_ave
               fn_flux(4, i, j) = theta*fn_flux(4, i, j) + (1.-theta)*e_ave
            else
            end if

			end if

         END DO
      END DO

      DO i = 0, nx + 1
         do j = 0, ny + 2
            den = u(1, i, j)
            xmt = u(2, i, j)
            ymt = u(3, i, j)
            eng = u(4, i, j)
            t0 = 1./den
            vex = xmt*t0
            vey = ymt*t0
            pre = gamma1*(eng - 0.5*(xmt*vex + ymt*vey))
            w(j) = sqrt(abs(den))
            vx(j) = vex
            vy(j) = vey
            h(j) = (pre + eng)*t0
         end do
         DO j = 0, ny + 1
            !********************************
            ! Roe Average Matrix for g'(u)
            !********************************

            t0 = w(j)/(w(j) + w(j + 1) + 1.e-15)
            t1 = 1.-t0
            vxm = t0*vx(j) + t1*vx(j + 1)
            vym = t0*vy(j) + t1*vy(j + 1)
            hm = t0*h(j) + t1*h(j + 1)
            qm = 0.5*(vxm*vxm + vym*vym)
            cm = sqrt(abs(gamma1*(hm - qm))) + 1.e-15
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

!          evr=0.
!          evl=0.
!          do ig=1,4
!             evr(ig,ig)=1.
!             evl(ig,ig)=1.
!          enddo

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
            !*****************************
            ! positivity limiter
            !*****************************
            !if (gp(1,i,j)<0. .or. gp(4,i,j)-gp(2,i,j)**2/gp(1,i,j)/2.-gp(3,i,j)**2/gp(1,i,j)/2.<0.) then
            !else
            !endif

            !if (gn(1,i,j)<0. .or. gn(4,i,j)-gn(2,i,j)**2/gn(1,i,j)/2.-gn(3,i,j)**2/gn(1,i,j)/2.<0.) then
            !else
            !endif


			if (use_limiter == 1) then
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
               p_l = gamma1*(e_l - m_l**2/d_l/2.-n_l**2/d_l/2.)
               if (p_l < 0.) theta = min(theta, abs(p_ave - 0.)/abs(p_ave - p_l))
            end do

            theta = max(theta, 0.)
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
               p_l = gamma1*(e_l - m_l**2/d_l/2.-n_l**2/d_l/2.)
               if (p_l < 0.) theta = min(theta, abs(p_ave - 0.)/abs(p_ave - p_l))
            end do
            theta = max(theta, 0.)
            if (theta < 1d0) then
               pre_count = pre_count + 1
               gn_flux(1, i, j) = theta*gn_flux(1, i, j) + (1.-theta)*d_ave
               gn_flux(2, i, j) = theta*gn_flux(2, i, j) + (1.-theta)*m_ave
               gn_flux(3, i, j) = theta*gn_flux(3, i, j) + (1.-theta)*n_ave
               gn_flux(4, i, j) = theta*gn_flux(4, i, j) + (1.-theta)*e_ave
            else
            end if

	     end if
         END DO
      END DO

   END SUBROUTINE wenofd
   !******************************
   SUBROUTINE scheme
      IMPLICIT NONE
      RHS = 0
      DO i = 1, nx
      DO j = 1, ny
         RHS(:, i, j) = -alphamax*((fp_flux(:, i, j) - fn_flux(:, i + 1, j)) - (fp_flux(:, i - 1, j) - fn_flux(:, i, j)))/dx&
         &        - alphamax*((gp_flux(:, i, j) - gn_flux(:, i, j + 1)) - (gp_flux(:, i, j - 1) - gn_flux(:, i, j)))/dy

      END DO
      END DO

   END SUBROUTINE scheme
   !**********************
   SUBROUTINE findalpha
      IMPLICIT NONE
      REAL(8)::den, pres, sound, v_x,v_y

      !alphamax = 0.
      !do i = -2, nx + 3
      !do j = -2, ny + 3
      !   den = u(1, i, j)
      !   pres = gamma1*(u(4, i, j) - u(2, i, j)**2/u(1, i, j)/2.-u(3, i, j)**2/u(1, i, j)/2.)
      !   sound = sqrt(gamma*pres/den)
      !   alphamax = max(alphamax, abs(u(2, i, j)/u(1, i, j)) + sound, abs(u(3, i, j)/u(1, i, j)) + sound)
      !end do
      !end do
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
         

   END SUBROUTINE findalpha
   !**************************
   SUBROUTINE rk
      IMPLICIT NONE
      REAL(8):: cfl, alphatemp,dt_cfl_x, dt_cfl_y, dt_cfl
      rktag = 1
      CALL wenofd

      cfl = 1./4.
      !dt = min(dx**(5./3.), dx*cfl/alphamax)
      !dt = dx**(5./3.)*cfl/alphamax

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
	   if (use_limiter == 2) then
         if (postag == 1) then
         call sweep_lim
         ! call scalar_sweep
         end if
      end if

	  CALL bc
      !if (postag==1) goto 110

      rktag = 2
      CALL wenofd

      CALL scheme
      u = 0.75d0*u0 + 0.25d0*(u + dt*RHS)
      call checkpositivity
      ! if (postag==1) goto 110

	   if (use_limiter == 2) then
         if (postag == 1) then
         call sweep_lim
         ! call scalar_sweep
         end if
      end if

      CALL bc

      rktag = 3
      CALL wenofd

      CALL scheme
      u = (u0 + 2d0*(u + dt*RHS))/3d0
      call checkpositivity
      !if (postag==1) goto 110
	   if (use_limiter == 2) then
         if (postag == 1) then
         call sweep_lim
         ! call scalar_sweep
         end if
      end if
      CALL bc

   END SUBROUTINE rk
   !*************************
   SUBROUTINE checkpositivity
      IMPLICIT NONE

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
   END SUBROUTINE checkpositivity
   !*************************
   subroutine bc
      implicit none
      !do j = -2, ny + 3
      !do i = 1, 3
      !   u(:, 1 - i, j) = u(:, nx + 1 - i, j)
      !   u(:, nx + i, j) = u(:, i, j)
      !end do
      !end do

      !do j = 1, 3
      !do i = -2, nx + 3
      !   u(:, i, 1 - j) = u(:, i, ny + 1 - j)
      !   u(:, i, ny + j) = u(:, i, j)
      !end do
      !end do

       do j = -2, ny + 3
            if (y(j) <= 0.05) then
               do i = 1, 3
                  u(1, 1 - i, j) = 5.
                  u(2, 1 - i, j) = 5.*ujet
                  u(3, 1 - i, j) = 0.
                  u(4, 1 - i, j) = 0.4127/gamma1 + 0.5*u(2, 1 - i, j)**2/u(1, 1 - i, j)
               end do
            else
               do i = 1, 3
                  u(1, 1 - i, j) = 0.5
                  u(2, 1 - i, j) = 0.
                  u(3, 1 - i, j) = 0.
                  u(4, 1 - i, j) = 0.4127/gamma1
               end do
            end if

         end do

      do j = 1, 3
         do i = -2, nx + 3
            u(:, i, 1 - j) = u(:, i, j)
            u(3, i, 1 - j) = -u(3, i, j)
         end do
      end do

      !if (x(nx + 1) > xm) then
         do ig = 1, 3
            u(:, nx + ig, :) = u(:, nx, :)
         end do

      !else
      !end if
      !if (y(ny + 1) > ym) then
         do ig = 1, 3
            u(:, :, ny + ig) = u(:, :, ny)
         end do
     ! else
     !end if
      return
   end subroutine bc
   !**************************
   subroutine initdata
      implicit none
      REAL(8):: xleft, xright, yleft, yright, ue(4)
     ! real(8):: rho_L, v_L, p_L, rho_R, v_R, p_R
      xleft = 0.
      xright = 1.
      dx = (xright - xleft)/nx

      yleft = 0.
      yright = 0.25
      dy = (yright - yleft)/ny

      do i = -2, nx + 3
         x(i) = xleft + i*dx - dx/2.
      end do
      do i = -2, ny + 3
         y(i) = yleft + i*dy - dy/2.
      end do

      u = 0.
      u(1, :, :) = 0.5
      u(4, :, :) = 0.4127/gamma1

      call bc
   end subroutine initdata
   !******************************
   FUNCTION exact(x, y, t)
      IMPLICIT NONE
      REAL(8), INTENT(IN)::x, y, t
      REAL(8), DIMENSION(1:4):: exact
      REAL(8):: mag, r
      

      r = sqrt((x - 5.-t)**2 + (y - 5.-t)**2)
      mag = 10.05
      exact(1) = (1 - gamma1*mag**2/8./gamma/pi**2*exp(1.-r**2))**(1./gamma1)
      exact(2) = exact(1)*(1.-mag/2./pi*exp((1.-r**2)/2.)*(y - 5.))
      exact(3) = exact(1)*(1.+mag/2./pi*exp((1.-r**2)/2.)*(x - 5.))
      exact(4) = exact(1)**gamma/gamma1 + 0.5*exact(2)**2/exact(1) + 0.5*exact(3)**2/exact(1)

!       exact(1)=1.+0.99*sin(2*pi*(x+y-2*t))
!       exact(2)=exact(1)
!       exact(3)=exact(1)
!       exact(4)=1./gamma1+exact(1)

   END FUNCTION exact
   !******************************
    SUBROUTINE output
      IMPLICIT NONE
      real(8):: lrho, lm, ln, le, pres
      

      open (73, file='den.dat')
      open(74, file='pre.dat')
      !open(75, file='vel1.dat')
      !open(76, file='vel2.dat')
      open(75, file='m1.dat')      
  !    open(74, file='pre.dat')

      !write (73, *) 'VARIABLES = "X", "y", "Numerical","Exact", "Error" '
     ! write (73, *) 'zone i=', nx, ', j=', ny
      write (73, *) 'VARIABLES = "X", "Y", "DENSITY"'

         write (73, *) 'zone i=', total_x, ', j=', total_y*2
         write (74, *) 'zone i=', total_x, ', j=', total_y*2
         write (75, *) 'zone i=', total_x, ', j=', total_y*2

      do j = 1, total_y
        do i = 1, total_x


           lrho = u(1, i, total_y - j + 1)


           lm = u(2, i, total_y - j + 1)
           ln = u(3, i, total_y - j + 1)
           le = u(4, i, total_y - j + 1)
           pres = gamma1*(le - 0.5*lm**2/lrho - 0.5*ln**2/lrho)


           write (73, *) dx*(i - 0.5), dy*(j - 0.5) - 0.25, log(abs(lrho))
           write (74, *) dx*(i - 0.5), dy*(j - 0.5) - 0.25, log(abs(pres))
           write (75, *) dx*(i - 0.5), dy*(j - 0.5) - 0.25, log(abs(pres/lrho))
        end do
        end do


        do j = 1, total_y
        do i = 1, total_x


           lrho = u(1, i, j)


           lm = u(2, i, j)
           ln = u(3, i, j)
           le = u(4, i, j)
           pres = gamma1*(le - 0.5*lm**2/lrho - 0.5*ln**2/lrho)


           write (73, *) dx*(i - 0.5), dy*(j - 0.5), log(abs(lrho)) !Shift to cell center (finite volume)
           write (74, *) dx*(i - 0.5), dy*(j - 0.5), log(abs(pres))
           write (75, *) dx*(i - 0.5), dy*(j - 0.5), log(abs(pres/lrho))
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
   END SUBROUTINE output

 !Sweeping limiter routines added by Chloe Griffin

   !*****************************
   FUNCTION nod_lim(v_1,v_2,v_bar)
      IMPLICIT NONE
      REAL(8), DIMENSION(4), INTENT(IN):: v_1,v_2,v_bar
      REAL(8):: p_v_1, p_v_2, p_v_bar, t,t_1, delta, d_1_2, d_1_bar
      REAL(8), DIMENSION(4):: nod_lim

    delta = 1d-6;


    P_v_1 = (gamma1)*(v_1(4) - 0.5*(v_1(2)*v_1(2) + v_1(3)*v_1(3))/v_1(1));
    P_v_2 = (gamma1)*(v_2(4) - 0.5*(v_2(2)*v_2(2) + v_2(3)*v_2(3))/v_2(1));
   

    !if (p_v_2 > 0 .AND. v_2(1)>0) then
    if (p_v_2 > eps) then
        t = ((P_v_1-eps)/(P_v_1 - P_v_2))
        t_pos_count = t_pos_count + 1 

    else
        P_v_bar = (gamma1)*(v_bar(4) - 0.5*(v_bar(2)*v_bar(2) + v_bar(3)*v_bar(3))/v_bar(1));
       ! t_1 = P_v_1/(P_v_1 - P_v_bar+1e-15);
        t_1 = (P_v_1-eps)/(P_v_1 - P_v_bar);
        d_1_2 = sqrt((v_2(4) - v_1(4))**2 + (v_2(3) - v_1(3))**2 + (v_2(2) - v_1(2))**2 + (v_2(1) - v_1(1))**2);
        d_1_bar = sqrt((v_bar(4) - v_1(4))**2 + (v_bar(3) - v_1(3))**2 + (v_bar(2) - v_1(2))**2 + (v_bar(1) - v_1(1))**2);
        
        !t = min(t_1*d_1_bar/d_1_2,1-delta);
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
   SUBROUTINE sweep_lim
   IMPLICIT NONE
   REAL(8), dimension(4)::  v_i, v_p, u_bar,  v_i_t, v_p_t, d
   !REAL(8), dimension(1:total_x, 1:total_y):: rho_t, m1_t, m2_t, E_t, P_check
   REAL(8), ALLOCATABLE :: rho_t(:,:), m1_t(:,:), m2_t(:,:), E_t(:,:)!P_check(:,:)
   REAL(8):: rho_avg, m1_avg, m2_avg, E_avg, w_p, pres_v_i , w_i, rho_min, p_min
   INTEGER:: q
   
   ALLOCATE(rho_t(1:total_x, 1:total_y))
   ALLOCATE(m1_t(1:total_x, 1:total_y))
   ALLOCATE(m2_t(1:total_x, 1:total_y))
   ALLOCATE(E_t(1:total_x, 1:total_y))
  ! ALLOCATE(P_check(1:total_x, 1:total_y))
if (use_limiter == 2) then   
   
   call scalar_sweep

   DO j = 1, total_y
      DO i = 1, total_x
            rho_t(i,j) = u(1,i,j);
            m1_t(i,j)  = u(2,i,j);
            m2_t(i,j)  = u(3,i,j);
            E_t(i,j)   = u(4,i,j);
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


   sweep_count = sweep_count + 1
   t_pos_count = 0
   t_neg_count = 0

    DO q = 1, swp_num
      if (mod(q,2) == 1) then  !Alternate sweep 1 and sweep 2
         !call scalar_sweep

        
        
            !Forward Sweep
            DO i = 1, (total_x*total_y-1)
      
                  v_i(1) = u(1,sweepI_i(i),sweepI_j(i));
                  v_i(2) = u(2,sweepI_i(i),sweepI_j(i));
                  v_i(3) = u(3,sweepI_i(i),sweepI_j(i));
                  v_i(4) = u(4,sweepI_i(i),sweepI_j(i));
      
                  pres_v_i = gamma1*(v_i(4) - 0.5*(v_i(2)*v_i(2) + v_i(3)*v_i(3))/(v_i(1)))
                  w_i = weights(sweepI_i(i),sweepI_j(i));
      
                  if (pres_v_i <= eps .OR. v_i(1) <= eps) then
      
                        v_p(1) = u(1,sweepI_i(i+1),sweepI_j(i+1));
                        v_p(2) = u(2,sweepI_i(i+1),sweepI_j(i+1));
                        v_p(3) = u(3,sweepI_i(i+1),sweepI_j(i+1));
                        v_p(4) = u(4,sweepI_i(i+1),sweepI_j(i+1));
                        
                        w_p = weights(sweepI_i(i+1),sweepI_j(i+1));
                        v_i_t = nod_lim(v_i,v_p,u_bar);
      
                           d  = v_i - v_i_t;
                           v_p_t =  v_p  + (w_i/w_p)*(d);
                           
                           u(1,sweepI_i(i),sweepI_j(i)) = v_i_t(1);
                           u(2,sweepI_i(i),sweepI_j(i)) = v_i_t(2);
                           u(3,sweepI_i(i),sweepI_j(i)) = v_i_t(3);
                           u(4,sweepI_i(i),sweepI_j(i)) = v_i_t(4);
      
                           u(1,sweepI_i(i+1),sweepI_j(i+1)) = v_p_t(1);
                           u(2,sweepI_i(i+1),sweepI_j(i+1)) = v_p_t(2);
                           u(3,sweepI_i(i+1),sweepI_j(i+1)) = v_p_t(3);
                           u(4,sweepI_i(i+1),sweepI_j(i+1)) = v_p_t(4);


                           ! if (i == 66 .AND. y_iter == 64) then
                           !end if
                  end if 
            end do
         
            !Backward Sweep
            DO i = total_x*total_y, 2, -1
            
                  v_i(1) = u(1,sweepI_i(i),sweepI_j(i));
                  v_i(2) = u(2,sweepI_i(i),sweepI_j(i));
                  v_i(3) = u(3,sweepI_i(i),sweepI_j(i));
                  v_i(4) = u(4,sweepI_i(i),sweepI_j(i));
      
                  pres_v_i =  gamma1*(v_i(4) - 0.5*(v_i(2)*v_i(2) + v_i(3)*v_i(3))/v_i(1));


                  w_i = weights(sweepI_i(i),sweepI_j(i));
      
                  if (pres_v_i <= eps  .OR. v_i(1) <= eps) then
      
                     v_p(1) = u(1,sweepI_i(i-1),sweepI_j(i-1));
                     v_p(2) = u(2,sweepI_i(i-1),sweepI_j(i-1));
                     v_p(3) = u(3,sweepI_i(i-1),sweepI_j(i-1));
                     v_p(4) = u(4,sweepI_i(i-1),sweepI_j(i-1));
                     
                     w_p = weights(sweepI_i(i-1),sweepI_j(i-1));
                     v_i_t = nod_lim(v_i,v_p,u_bar);
      
                           d  = v_i - v_i_t;
      
                           v_p_t =  v_p  + (w_i/ w_p)*(d);
                           
                           u(1,sweepI_i(i),sweepI_j(i)) = v_i_t(1);
                           u(2,sweepI_i(i),sweepI_j(i)) = v_i_t(2);
                           u(3,sweepI_i(i),sweepI_j(i)) = v_i_t(3);
                           u(4,sweepI_i(i),sweepI_j(i)) = v_i_t(4);
      
                           u(1,sweepI_i(i-1),sweepI_j(i-1)) = v_p_t(1);
                           u(2,sweepI_i(i-1),sweepI_j(i-1)) = v_p_t(2);
                           u(3,sweepI_i(i-1),sweepI_j(i-1)) = v_p_t(3);
                           u(4,sweepI_i(i-1),sweepI_j(i-1)) = v_p_t(4);
      
                  end if
            end do


   else

      
          !Forward Sweep
           DO i = 1, (total_x*total_y-1)
    
                 v_i(1) = u(1,sweepII_i(i),sweepII_j(i));
                 v_i(2) = u(2,sweepII_i(i),sweepII_j(i));
                 v_i(3) = u(3,sweepII_i(i),sweepII_j(i));
                 v_i(4) = u(4,sweepII_i(i),sweepII_j(i));
    
                 pres_v_i = gamma1*(v_i(4) - 0.5*(v_i(2)*v_i(2) + v_i(3)*v_i(3))/(v_i(1)))
                 w_i = weights(sweepII_i(i),sweepII_j(i));
    
                 if (pres_v_i <= eps .OR. v_i(1) <= eps) then
    
                       v_p(1) = u(1,sweepII_i(i+1),sweepII_j(i+1));
                       v_p(2) = u(2,sweepII_i(i+1),sweepII_j(i+1));
                       v_p(3) = u(3,sweepII_i(i+1),sweepII_j(i+1));
                       v_p(4) = u(4,sweepII_i(i+1),sweepII_j(i+1));
                      
                       w_p = weights(sweepII_i(i+1),sweepII_j(i+1));
                       v_i_t = nod_lim(v_i,v_p,u_bar);
    
                          d  = v_i - v_i_t;
                          v_p_t =  v_p  + (w_i/w_p)*(d);
                         
                          u(1,sweepII_i(i),sweepII_j(i)) = v_i_t(1);
                          u(2,sweepII_i(i),sweepII_j(i)) = v_i_t(2);
                          u(3,sweepII_i(i),sweepII_j(i)) = v_i_t(3);
                          u(4,sweepII_i(i),sweepII_j(i)) = v_i_t(4);
    
                          u(1,sweepII_i(i+1),sweepII_j(i+1)) = v_p_t(1);
                          u(2,sweepII_i(i+1),sweepII_j(i+1)) = v_p_t(2);
                          u(3,sweepII_i(i+1),sweepII_j(i+1)) = v_p_t(3);
                          u(4,sweepII_i(i+1),sweepII_j(i+1)) = v_p_t(4);




                          ! if (i == 66 .AND. y_iter == 64) then
                          !end if
                 end if
           end do
       
           !Backward Sweep
           DO i = total_x*total_y, 2, -1
          
                 v_i(1) = u(1,sweepII_i(i),sweepII_j(i));
                 v_i(2) = u(2,sweepII_i(i),sweepII_j(i));
                 v_i(3) = u(3,sweepII_i(i),sweepII_j(i));
                 v_i(4) = u(4,sweepII_i(i),sweepII_j(i));
    
                 pres_v_i =  gamma1*(v_i(4) - 0.5*(v_i(2)*v_i(2) + v_i(3)*v_i(3))/v_i(1));




                 w_i = weights(sweepII_i(i),sweepII_j(i));
    
                 if (pres_v_i <= eps  .OR. v_i(1) <= eps) then
    
                    v_p(1) = u(1,sweepII_i(i-1),sweepII_j(i-1));
                    v_p(2) = u(2,sweepII_i(i-1),sweepII_j(i-1));
                    v_p(3) = u(3,sweepII_i(i-1),sweepII_j(i-1));
                    v_p(4) = u(4,sweepII_i(i-1),sweepII_j(i-1));
                   
                    w_p = weights(sweepII_i(i-1),sweepII_j(i-1));
                    v_i_t = nod_lim(v_i,v_p,u_bar);
    
                          d  = v_i - v_i_t;
    
                          v_p_t =  v_p  + (w_i/ w_p)*(d);
                         
                          u(1,sweepII_i(i),sweepII_j(i)) = v_i_t(1);
                          u(2,sweepII_i(i),sweepII_j(i)) = v_i_t(2);
                          u(3,sweepII_i(i),sweepII_j(i)) = v_i_t(3);
                          u(4,sweepII_i(i),sweepII_j(i)) = v_i_t(4);
    
                          u(1,sweepII_i(i-1),sweepII_j(i-1)) = v_p_t(1);
                          u(2,sweepII_i(i-1),sweepII_j(i-1)) = v_p_t(2);
                          u(3,sweepII_i(i-1),sweepII_j(i-1)) = v_p_t(3);
                          u(4,sweepII_i(i-1),sweepII_j(i-1)) = v_p_t(4);
    
                 end if
           end do
         end if 
         rho_min = 1d20
         P_min = 1d20
         DO j = 1, total_y
            DO i = 1, total_x
               !P_check(j) = gamma1*(u(4,x_iter,j) - 0.5*(u(2,x_iter,j)**2 + u(3,x_iter,j)**2)/u(1,x_iter,j));
               P_min = min(gamma1*(u(4,i,j) - 0.5*(u(2,i,j)*u(2,i,j) + u(3,i,j)*u(3,i,j))/u(1,i,j)),P_min);
               rho_min = min(u(1,i,j), rho_min);
            end do
         end do
   
         pre_count = pre_count + 1  
         

         if (P_min >= eps .AND. rho_min>=eps) then !(P_min > 0 .AND. rho_min>0) then
            sweep_track(sweep_count) = q
            neg_track(sweep_count) = t_neg_count
            pos_track(sweep_count) = t_pos_count

            exit
         end if

         if (q == swp_num) then
            print *, "SWEEP FAIL: swp_num =",swp_num
         end if

      end do
   end if 
  
      
   postag = 0

      END SUBROUTINE  sweep_lim

!*****************************
      SUBROUTINE scalar_sweep
         IMPLICIT NONE
         REAL(8)::    v_i, v_ip1, v_im1
         REAL(8), ALLOCATABLE :: rho_t(:,:)
         !REAL(8), dimension(1:total_x,1:total_y):: rho_t 
         REAL(8)::  w_i, w_ip1, w_im1, delta, rho_min, rho_max

       ALLOCATE(rho_t(1:total_x, 1:total_y))  
               
      if (use_limiter == 2) then
      rho_min = 1d20
      DO j = 1, total_y
         DO i = 1, total_x
             rho_min = min(u(1,i,j), rho_min);
         end do
      end do

      delta = 1d-12;
      if (rho_min < delta) then
      DO j = 1, total_y
         DO i = 1, total_x
            rho_t(i,j) = u(1,i,j);
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
            u(1,i,j) = rho_t(i,j);
         end do
      end do

      den_count = den_count + 1;
   end if 

      rho_min = 1d20
      DO j = 1, total_y
         DO i = 1, total_x
            rho_min = min(u(1,i,j), rho_min);
         end do
      end do

      rho_max = 0
      DO j = 1, total_y
         DO i = 1, total_x
            rho_max = max(u(1,i,j), rho_max);
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



END PROGRAM FDWENO4EULER1D
