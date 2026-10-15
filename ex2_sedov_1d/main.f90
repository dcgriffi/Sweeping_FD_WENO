PROGRAM FDWENO4EULER1D
   USE weno5
   IMPLICIT NONE
   INTEGER, PARAMETER:: n = 800, ierror = 0, USE_LIMITER = 2, RK_choice = 1, swp_num = 50
   REAL, PARAMETER::  nSteps = 800, gamma=1.4, gamma1=gamma-1., eps=1d-13, time_terminal=0.001
 !  REAL, PARAMETER:: factor = 3, cfl=0.4*factor, max_iter = 10000
    REAL, PARAMETER:: cfl=1.2, max_iter = 10000

!factor = 3.2
   REAL, DIMENSION(3, -2:n + 3):: u, u0, RHS, flus, hp_flux, hn_flux, hp, hn
   REAL, DIMENSION(-2:n + 3):: x, weights
   integer:: rktag, kcount
   real:: dx, dt, time, alphamax
   integer:: i, j, q, k, kk, m, den_count = 0, pre_count = 0, postag

    !Adding sweep data
   INTEGER, PARAMETER :: total_max = 5000
   INTEGER, DIMENSION(1:total_max) :: sweep_track = 0
   INTEGER :: sweep_count = 0, total_sweep, max_sweep
   real(8):: avg_sweep

   INTEGER, DIMENSION(1:total_max) :: neg_track = 0
   INTEGER, DIMENSION(1:total_max) :: pos_track = 0
   INTEGER :: t_neg_count, t_pos_count, total_neg, total_pos, max_neg, max_pos
   real(8):: avg_pos, avg_neg


   call initdata

   DO i = -2, n + 3
      weights(i) = 1/(nSteps + 6)
   END DO


   time = 0d0
   kcount = 0

   DO WHILE (time <= time_terminal - 1d-10)
      if (RK_choice == 1) then
         CALL rk3
      else
         call rk
      end if
      time = time + dt
      kcount = kcount + 1
      if (kcount/50*50 .eq. kcount) write (6, *) kcount, time/time_terminal*100d0, '%'

      if (kcount .eq. max_iter) then
         print *, "min_rho:", minval(u(1, :))
         CALL output
         stop "Finished to max_iters"
      end if
   END DO
   CALL output

contains
   FUNCTION fflux(rho, mon_x, e)
      IMPLICIT NONE
      REAL, INTENT(IN):: rho, mon_x, e
      REAL:: vx, pres
      REAL, DIMENSION(1:3)::fflux

      if (rho > 1d-8) vx = mon_x/rho
      if (rho <= 1d-8) vx = 0.

      pres = gamma1*(e - 0.5*rho*vx**2)

      fflux(1) = mon_x
      fflux(2) = rho*vx**2 + pres
      fflux(3) = (e + pres)*vx

      RETURN
   END FUNCTION fflux
   !**************************
   SUBROUTINE wenofd
      IMPLICIT NONE
      REAL:: evr(3, 3), evl(3, 3), cu_p(3, -2:2), cu_n(3, -2:2)
      REAL:: den, xmt, eng, t0, t1, t2, pre, vex, cvel, vxm, hm, qm, cm, rcm, b1, b2
      REAL:: ul(3), ur(3), d_ave, m_ave, e_ave, p_ave, theta, point(3, 2), d_l, m_l, e_l, p_l
      call findalpha

      do i = -2, n + 3
         hp(:, i) = 0.5*(u(:, i) + fflux(u(1, i), u(2, i), u(3, i))/alphamax)
         hn(:, i) = 0.5*(u(:, i) - fflux(u(1, i), u(2, i), u(3, i))/alphamax)
      end do

      call checkpositivity

      DO m = 0, n + 1
!          !********************************
!          ! Roe Average Matrix
!          !********************************
!          t0 =sqrt(u(1,m))/(sqrt(u(1,m))+sqrt(u(1,m+1)))
!          t1 = 1.-t0
!          vxm = t0*u(2,m)/u(1,m)+t1*u(2,m+1)/u(1,m+1)
!          hm = t0*(gamma1*(u(3,m)-0.5*u(2,m)**2/u(1,m))+u(3,m))/u(1,m) &
!             &+t1*(gamma1*(u(3,m+1)-0.5*u(2,m+1)**2/u(1,m+1))+u(3,m+1))/u(1,m+1)
!          qm = 0.5*vxm*vxm
!          cm = sqrt( gamma1*(hm-qm))
!          t0 = vxm*cm

         t0 = 0.5
         t1 = 0.5
         vxm = u(2, m)/u(1, m)
         hm = (gamma1*(u(3, m) - 0.5*u(2, m)**2/u(1, m)) + u(3, m))/u(1, m)
         qm = 0.5*vxm*vxm
         cm = sqrt(gamma1*(hm - qm))
         t0 = vxm*cm

         evr(1, 1) = 1.0
         evr(1, 2) = 1.0
         evr(1, 3) = 1.0
         evr(2, 1) = vxm - cm
         evr(2, 2) = vxm
         evr(2, 3) = vxm + cm
         evr(3, 1) = hm - t0
         evr(3, 2) = qm
         evr(3, 3) = hm + t0
         rcm = 1./cm
         b1 = gamma1*rcm*rcm
         b2 = qm*b1
         t0 = vxm*rcm
         t1 = b1*vxm
         t2 = 0.5*b1
         evl(1, 1) = 0.5*(b2 + t0)
         evl(1, 2) = -0.5*(t1 + rcm)
         evl(1, 3) = t2
         evl(2, 1) = 1.-b2
         evl(2, 2) = t1
         evl(2, 3) = -b1
         evl(3, 1) = 0.5*(b2 - t0)
         evl(3, 2) = -0.5*(t1 - rcm)
         evl(3, 3) = t2

         !******************************************
         ! projection into characteristic fields
         !******************************************
         DO j = -2, 2
         DO k = 1, 3
            cu_p(k, j) = evl(k, 1)*hp(1, m + j) + evl(k, 2)*hp(2, m + j) + evl(k, 3)*hp(3, m + j)
            cu_n(k, j) = evl(k, 1)*hn(1, m + j) + evl(k, 2)*hn(2, m + j) + evl(k, 3)*hn(3, m + j)
         END DO
         END DO
         !**************************
         ! weno reconstruction
         !**************************
         DO kk = 1, 3
            call weno_right(cu_p(kk, :), ur(kk))
            call weno_left(cu_n(kk, :), ul(kk))
         END DO
         !******************************************
         ! projection back into conserved variables
         !******************************************
         DO k = 1, 3
            hp_flux(k, m) = evr(k, 1)*ur(1) + evr(k, 2)*ur(2) + evr(k, 3)*ur(3)
            hn_flux(k, m) = evr(k, 1)*ul(1) + evr(k, 2)*ul(2) + evr(k, 3)*ul(3)
         END DO

         if (USE_LIMITER == 1)  then


            !if (postag == 1) then
               !*****************************
               ! positivity limiter
               !*****************************
               theta = 1d0
               d_ave = hp(1, m)
               m_ave = hp(2, m)
               e_ave = hp(3, m)
               p_ave = gamma1*(e_ave - 0.5*m_ave**2/d_ave)
               theta = 1.
               point(:, 1) = hp_flux(:, m)
               point(:, 2) = (hp(:, m) - hp_flux(:, m)/12.)/(11./12.)

               do i = 1, 2
                  d_l = point(1, i)
                  if (d_l < 0.) theta = min(theta, abs(d_ave - 0.)/abs(d_ave - d_l))
               end do
               if (theta < 1.) then
                  hp_flux(1, m) = theta*hp_flux(1, m) + (1.-theta)*d_ave
                  den_count = den_count + 1
               else
               end if

               theta = 1.
               point(1, 1) = hp_flux(1, m)
               point(1, 2) = (hp(1, m) - hp_flux(1, m)/12.)/(11./12.)

               do i = 1, 2
                  d_l = point(1, i)
                  m_l = point(2, i)
                  e_l = point(3, i)
                  p_l = gamma1*(e_l - m_l**2/d_l/2.)
                  if (p_l < 0.) theta = min(theta, abs(p_ave - 0.)/abs(p_ave - p_l))
               end do

               if (theta < 1d0) then
                  pre_count = pre_count + 1
                  hp_flux(1, m) = theta*hp_flux(1, m) + (1.-theta)*d_ave
                  hp_flux(2, m) = theta*hp_flux(2, m) + (1.-theta)*m_ave
                  hp_flux(3, m) = theta*hp_flux(3, m) + (1.-theta)*e_ave

               else
               end if

               theta = 1d0
               d_ave = hn(1, m)
               m_ave = hn(2, m)
               e_ave = hn(3, m)
               p_ave = gamma1*(e_ave - 0.5*m_ave**2/d_ave)
               theta = 1.
               point(:, 1) = hn_flux(:, m)
               point(:, 2) = (hn(:, m) - hn_flux(:, m)/12.)/(11./12.)

               do i = 1, 2
                  d_l = point(1, i)
                  if (d_l < 0.) theta = min(theta, abs(d_ave - 0.)/abs(d_ave - d_l))
               end do
               if (theta < 1.) then
                  hn_flux(1, m) = theta*hn_flux(1, m) + (1.-theta)*d_ave
                  den_count = den_count + 1
               else
               end if

               theta = 1.
               point(1, 1) = hn_flux(1, m)
               point(1, 2) = (hn(1, m) - hn_flux(1, m)/12.)/(11./12.)

               do i = 1, 2
                  d_l = point(1, i)
                  m_l = point(2, i)
                  e_l = point(3, i)
                  p_l = gamma1*(e_l - m_l**2/d_l/2.)
                  if (p_l < 0.) theta = min(theta, abs(p_ave - 0.)/abs(p_ave - p_l))
               end do

               if (theta < 1d0) then
                  pre_count = pre_count + 1
                  hn_flux(1, m) = theta*hn_flux(1, m) + (1.-theta)*d_ave
                  hn_flux(2, m) = theta*hn_flux(2, m) + (1.-theta)*m_ave
                  hn_flux(3, m) = theta*hn_flux(3, m) + (1.-theta)*e_ave

               else
               end if
            !end if 
         end if

      END DO

   END SUBROUTINE wenofd
   !******************************
   SUBROUTINE scheme
      IMPLICIT NONE
!       do i=1,n
!          if(i==99) then
!
!          else
!          endif
!       enddo

      DO i = 1, n + 1
         flus(:, i) = alphamax*(hp_flux(:, i - 1) - hn_flux(:, i))
      END DO

      RHS = 0d0 !Important fix here (by Chloe Griffin)

      DO i = 1, n
         RHS(:, i) = -(flus(:, i + 1) - flus(:, i))/dx
      END DO

   END SUBROUTINE scheme
   !**********************
   SUBROUTINE findalpha
      IMPLICIT NONE
      REAL::den, pres, sound

      alphamax = 0.
      do i = -2, n + 3
         den = u(1, i)
         pres = gamma1*(u(3, i) - u(2, i)**2/u(1, i)/2.)
         sound = sqrt(gamma*pres/den)
         alphamax = max(alphamax, abs(u(2, i)/u(1, i)) + sound)
      end do

   END SUBROUTINE findalpha
   !**************************
   SUBROUTINE rk3
      IMPLICIT NONE
      REAL:: alphatemp
      rktag = 1
      CALL wenofd

      dt = dx*cfl/alphamax

      alphatemp = alphamax
      postag = 0
      u0 = u
      CALL scheme
      u = u + dt*RHS
      call checkpositivity

      if (postag == 1) then
         CALL sweep_lim
      end if

      rktag = 2
      CALL wenofd

      CALL scheme
      u = 0.75d0*u0 + 0.25d0*(u + dt*RHS)
      call checkpositivity
      
      if (postag == 1) then
         CALL sweep_lim
      end if

      rktag = 3
      CALL wenofd

      CALL scheme
      u = (u0 + 2d0*(u + dt*RHS))/3d0
      call checkpositivity
      
      if (postag == 1) then
         CALL sweep_lim
      end if

   END SUBROUTINE rk3

   SUBROUTINE rk
      IMPLICIT NONE
      REAL:: alphatemp
      rktag = 1
      CALL wenofd

      dt = dx*cfl/alphamax

      alphatemp = alphamax
      postag = 0
110   continue
      rktag = 1
      if (postag == 1) then
         u = u0
         dt = dt/2.
         alphamax = alphatemp
      else
      end if

      u0 = u
      CALL scheme
      u = u + dt*RHS
      call checkpositivity
      if (postag == 1) goto 110
      CALL bc

      rktag = 2
      CALL wenofd

      CALL scheme
      u = 0.75d0*u0 + 0.25d0*(u + dt*RHS)
      call checkpositivity
      if (postag == 1) goto 110
      CALL bc

      rktag = 3
      CALL wenofd

      CALL scheme
      u = (u0 + 2d0*(u + dt*RHS))/3d0
      call checkpositivity
      if (postag == 1) goto 110
      CALL bc

   END SUBROUTINE rk
   !*************************
   SUBROUTINE checkpositivity
      IMPLICIT NONE

      postag = 0
      do m = 1, n
         if (u(1, m) <= eps) then
            postag = 1
         else
         end if

         !if (u(3, m) - u(2, m)**2/u(1, m)/2. < eps) then
         if (u(3, m) - u(2, m)**2/u(1, m)/2. <= eps) then
            postag = 1
         else
         end if
      end do
   END SUBROUTINE checkpositivity
   !*************************
   subroutine bc
      implicit none
!       u(:,0)=u(:,n)
!       u(:,-1)=u(:,n-1)
!       u(:,-2)=u(:,n-2)
!       u(:,n+1)=u(:,1)
!       u(:,n+2)=u(:,2)
!       u(:,n+3)=u(:,3)
      return
   end subroutine bc
   !**************************
   subroutine initdata
      implicit none
      real:: xleft, xright
      xleft = -2.
      xright = 2.
      dx = (xright - xleft)/n

      do i = -2, n + 3
         x(i) = xleft + i*dx - dx/2.
      end do

      do j = -2, n + 3
         u(:, j) = exact(x(j), 0.)
      end do

      u = 0d0
      do i = -2, n + 3
         u(1, i) = 1.
         u(3, i) = 1d-9/gamma1
      end do
      u(3, n/2 + 1) = 3200000.0/dx

      call bc
   end subroutine initdata
   !******************************
   FUNCTION exact(x, t)
      IMPLICIT NONE
      REAL, INTENT(IN)::x, t
      REAL, DIMENSION(1:3):: exact

      exact(1) = 1.+0.99999999*sin(2*pi*(x - t))
      exact(2) = exact(1)
      exact(3) = 1./gamma1 + 0.5*exact(1)

   END FUNCTION exact
   !******************************
   SUBROUTINE output
      IMPLICIT NONE
      REAL:: er1, eri, temp, eriold, er1old, o1, oi
      REAL:: uue(3), xrg

      er1 = 0.
      eri = 0.

      do i = 1, n
         uue = exact(x(i), time)
         er1 = er1 + abs(uue(1) - u(1, i))
         eri = max(eri, abs(uue(1) - u(1, i)))
      end do

      er1 = er1/n

      eriold = 1.
      er1old = 1.
      if (ierror .eq. 1) then
         read (7, 11) er1old, eriold
      end if

      o1 = log(er1old/er1)/log(2D0)
      oi = log(eriold/eri)/log(2D0)

      print *, 'L1 = ', ' Linfinity ='
      print *, er1, eri
      print *, 'L1order=', ' Liorder='
      print *, o1, oi

11    format(6e13.5)
      rewind (7)
      write (7, 11) er1, eri

      open (73, file='den.dat')
      open (74, file='pre.dat')
      open (75, file='vel.dat')

      do i = 1, n
         write (73, *) x(i), u(1, i)
         write (74, *) x(i), gamma1*(u(3, i) - 0.5*u(2, i)**2/u(1, i))
         write (75, *) x(i), u(2, i)/u(1, i)
      end do

      total_sweep = sum(sweep_track(1:sweep_count))
      avg_sweep = real(total_sweep,8)/real(sweep_count,8)
      max_sweep = maxval(sweep_track(1:sweep_count))

      print *, 'total sweep = ', total_sweep
      print *, 'avg sweep = ',   avg_sweep
      print *, 'max sweep = ',   max_sweep

      !Sweep count is number of times sweeping is triggered

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
      REAL, DIMENSION(3), INTENT(IN):: v_1,v_2,v_bar
      REAL:: p_v_1, p_v_2, p_v_bar, t,t_1, delta, d_1_2, d_1_bar
      REAL, DIMENSION(3)::nod_lim


    P_v_1 = (gamma1)*(v_1(3) - 0.5*(v_1(2)**2/v_1(1)));
    P_v_2 = (gamma1)*(v_2(3) - 0.5*(v_2(2)**2/v_2(1)));
   

    !if (p_v_2 > 0 .AND. v_2(1)>0) then
    if (p_v_2 > eps) then

        t = (P_v_1-eps)/(P_v_1 - P_v_2);
        t_pos_count = t_pos_count + 1
    else
        P_v_bar = (gamma1)*(v_bar(3) - 0.5*(v_bar(2)**2/v_bar(1)));


        t_1 = (P_v_1-eps)/(P_v_1 - P_v_bar);

        d_1_2 = sqrt((v_2(3) - v_1(3))**2 + (v_2(2) - v_1(2))**2 + (v_2(1) - v_1(1))**2);
        d_1_bar = sqrt((v_bar(3) - v_1(3))**2 + (v_bar(2) - v_1(2))**2 + (v_bar(1) - v_1(1))**2);
        
        t = min(t_1*d_1_bar/d_1_2,0.25);

        t_neg_count = t_neg_count + 1
     
    end if






    nod_lim(1) = (1 - t)*v_1(1) + t*v_2(1);
    nod_lim(2) = (1 - t)*v_1(2) + t*v_2(2);
    nod_lim(3) = (1 - t)*v_1(3) + t*v_2(3);



      RETURN
   END FUNCTION nod_lim


   !*****************************
   SUBROUTINE sweep_lim
   IMPLICIT NONE
   REAL, dimension(3)::  v_j, v_p, u_bar, v_j_t, v_p_t, d
   REAL, dimension(-2:n + 3):: rho_t, m_t, E_t, P_check
   REAL:: rho_avg, m_avg, E_avg, pres_v_j , w_p, w_j, rho_min, p_min
   
   if (use_limiter == 2) then

   call scalar_sweep

   DO i = -2,  n + 3
         rho_t(i) = u(1,i);
         m_t(i) = u(2,i);
         E_t(i) = u(3,i);
   END DO

      
   rho_avg = dot_product(rho_t,weights);
   m_avg = dot_product(m_t,weights); !m_t*w';
   E_avg = dot_product(E_t,weights); !E_t*w';
      
      
      u_bar(1) = rho_avg
      u_bar(2) =  m_avg 
      u_bar(3) = E_avg;


      sweep_count = sweep_count + 1
      t_pos_count = 0
      t_neg_count = 0
      DO q = 1, swp_num
              !Forward Sweep
              DO j = -2, n + 2
      
      
                  v_j(1) = u(1,j);
                  v_j(2) = u(2,j);
                  v_j(3) = u(3,j);
      
                  pres_v_j = gamma1*(v_j(3) - 0.5*(v_j(2)**2)/v_j(1));
      
      
              
                  
                  
                  w_j = weights(j); 
      
      
                  if (pres_v_j <= eps .OR. v_j(1) <= eps) then
      
                         
                         v_p(1) = u(1,j+1);
                         v_p(2) = u(2,j+1);
                         v_p(3) = u(3,j+1);
                        
                        
                        w_p = weights(j+1);
                        v_j_t = nod_lim(v_j,v_p,u_bar);
      
             
      
                            d  = v_j - v_j_t; 
                            v_p_t =  v_p  + (w_j/w_p)*(d); 


                          ! if (j == 401) then
                          ! end if 
                            
                            u(1,j) = v_j_t(1);
                            u(2,j) = v_j_t(2);
                            u(3,j) = v_j_t(3);
      
                            u(1,j+1) = v_p_t(1);
                            u(2,j+1) = v_p_t(2);
                            u(3,j+1) = v_p_t(3);
      
                  end if
               
      
      
               end do
      
      
       
           
          
              !Backward Sweep
              DO j = n + 3, -1, -1
             
      
                  v_j(1) = u(1,j);
                  v_j(2) = u(2,j);
                  v_j(3) = u(3,j);
      
                  pres_v_j = gamma1*(v_j(3) - 0.5*(v_j(2)**2)/v_j(1));
      
      
                  w_j = weights(j);
      
          
      
      
                  if (pres_v_j <= eps  .OR. v_j(1) <= eps) then
      
                     v_p(1) = u(1,j-1); 
                     v_p(2) = u(2,j-1); 
                     v_p(3) = u(3,j-1); 
                      
                      
                      w_p = weights(j-1); 
                      v_j_t = nod_lim(v_j,v_p,u_bar);
      
            
      
                           d  = v_j - v_j_t; 
      
                           v_p_t =  v_p  + (w_j/ w_p)*(d);

                         

                          ! if (j == 401) then
                          ! end if 
                            
                           u(1,j) = v_j_t(1);
                           u(2,j) = v_j_t(2);
                           u(3,j) = v_j_t(3);
      
                           u(1,j-1) = v_p_t(1);
                           u(2,j-1) = v_p_t(2);
                           u(3,j-1) = v_p_t(3);
       
      
                  end if
              end do
            
              rho_min = 1d20
              P_min = 1d20
              DO i = -2, n + 3
                  P_check(i) = gamma1*(u(3,i) - 0.5*(u(2,i)**2)/u(1,i));
                  P_min = min(gamma1*(u(3,i) - 0.5*(u(2,i)**2)/u(1,i)),P_min);
                  rho_min = min(u(1,i), rho_min);
              end do

      
              pre_count = pre_count + 1
              
              
              if (P_min > eps .AND. rho_min>eps) then !(P_min > 0 .AND. rho_min>0) then
               postag = 0
               sweep_track(sweep_count) = q
               neg_track(sweep_count) = t_neg_count
               pos_track(sweep_count) = t_pos_count
               exit
              end if


   
      
      end do
   end if
      

      END SUBROUTINE  sweep_lim

!*****************************
      SUBROUTINE scalar_sweep
         IMPLICIT NONE
         REAL::  v_j, v_jp1, v_jm1 
         REAL, dimension(-2:n + 3):: rho_t 
         REAL::  w_j, w_jp1, w_jm1, delta, rho_min, rho_max
               
    if (use_limiter == 2) then

         rho_min = 1d20
         DO i = -2, n + 3
             rho_min = min(u(1,i), rho_min);
         end do

      delta = 1d-13;

      DO j = -2, n + 3
         rho_t(j) = u(1,j);
      end do
      
      ! TO DO:
      !Check that array sizes match
      !Check that weights sum to one


      DO j = -2, n+2
         v_j = rho_t(j); v_jp1 = rho_t(j + 1); w_j = weights(j); w_jp1 = weights(j + 1);
         if (v_j < delta) then
            rho_t(j + 1) =  v_jp1  + (w_j/w_jp1)*(v_j - delta);
            rho_t(j) = delta;
            end if
      end do


      DO j = n + 3, -1, -1
         v_j = rho_t(j); v_jm1 = rho_t(j - 1); w_j = weights(j); w_jm1 = weights(j - 1);
         if (v_j < delta) then
            rho_t(j - 1) =  v_jm1  + (w_j/ w_jm1)*(v_j - delta);
            rho_t(j) = delta;
         end if
      end do

      DO  j = -2, n + 3
         u(1,j) = rho_t(j);
      end do

      den_count = den_count + 1;

      rho_min = 1d20
      DO i = -2, n + 3
          rho_min = min(u(1,i), rho_min);
      end do

      rho_max = 0
      DO i = -2, n + 3
          rho_max = max(u(1,i), rho_max);
      end do
      postag = 0

   end if
            
      
 END SUBROUTINE  scalar_sweep
      
      



END PROGRAM FDWENO4EULER1D
