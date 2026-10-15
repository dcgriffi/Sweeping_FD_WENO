      MODULE weno5
      IMPLICIT NONE
      REAL, PARAMETER, DIMENSION(-1:2, 0:2):: c=RESHAPE((/ &
      1.83333333333333333333333333333333,& 
      0.333333333333333333333333333333333,&
      -0.166666666666666666666666666666667,&
      0.333333333333333333333333333333333,&
      -1.16666666666666666666666666666667,&
      0.833333333333333333333333333333333,&
      0.833333333333333333333333333333333,&
      -1.16666666666666666666666666666667,&
      0.333333333333333333333333333333333,&
      -0.166666666666666666666666666666667,&
      0.333333333333333333333333333333333,&
      1.83333333333333333333333333333333/),(/4,3/))
      REAL, PARAMETER, DIMENSION(1:6):: xg=(/&
      -0.466234757101576013906150777246997,&
      -0.330604693233132256830699797509953,&
      -0.119309593041598454315250860840356,&
       0.119309593041598454315250860840356,&
       0.330604693233132256830699797509953,&
       0.466234757101576013906150777246997/), &
       wg=(/8.566224618958517252014807108636649E-0002,&
       0.180380786524069303784916756918858,&
       0.233956967286345523694935171994775,&
       0.233956967286345523694935171994775,&
       0.180380786524069303784916756918858,&
       8.566224618958517252014807108636649E-0002/)
       REAL, PARAMETER, DIMENSION(0:2, 0:2, 1:3):: cq=RESHAPE((/ &
        1.61428083526444586611022314330069,&
  	0.226982500643704177592296603322453,&
 	-0.160315833977037510925629936655787,&
 	-0.841263335908150043702519746623147,&
 	 0.933333333333333333333333333333333,&
  	0.707930002574816710369186413289813,&
  	0.226982500643704177592296603322453,&
 	-0.160315833977037510925629936655787,&
  	0.452385831402220800556443523365973,&
  	0.958333333333333333333333333333333,&
 	-4.166666666666666666666666666666666E-0002,&
 	-4.166666666666666666666666666666666E-0002,&
  	8.333333333333333333333333333333333E-0002,&
   	1.08333333333333333333333333333333,&
  	8.333333333333333333333333333333333E-0002,&
 	-4.166666666666666666666666666666666E-0002,&
 	-4.166666666666666666666666666666666E-0002,&
  	0.958333333333333333333333333333333,&
  	0.452385831402220800556443523365973,&
 	-0.160315833977037510925629936655787,&
  	0.226982500643704177592296603322453,&
  	0.707930002574816710369186413289813,&
  	0.933333333333333333333333333333333,&
 	-0.841263335908150043702519746623147,&
 	-0.160315833977037510925629936655787,&
  	0.226982500643704177592296603322453,&
   	1.61428083526444586611022314330069/),(/3,3,3/))
	!************************************************************************************
	! rq is the linear weights for the 3-point Gaussian quadrature points
	! rq(:,2) need to be dealt with the splitting technique to remove the negative weight
	!***********************************************************************************
	REAL, PARAMETER, DIMENSION(1:3, 1:3)::rq=RESHAPE((/ &
	  0.139888966110548358998525220728139       ,&
          0.615267175572519083969465648854962       ,&
          0.244843858316932557032009130416899       ,&
         -0.112500000000000000000000000000000       ,&
           1.22500000000000000000000000000000       ,&
         -0.112500000000000000000000000000000       ,&
          0.244843858316932557032009130416899       ,&
          0.615267175572519083969465648854962       ,&
          0.139888966110548358998525220728139/),(/3,3/))
	REAL, PARAMETER, DIMENSION(1:3):: rq_p=(/ &
	4.205607476635514018691588785046729E-0002,&
        0.915887850467289719626168224299065,&
        4.205607476635514018691588785046729E-0002/), &
        rq_n=(/0.134328358208955223880597014925373,&
        0.731343283582089552238805970149254,&
        0.134328358208955223880597014925373/)
        REAL, PARAMETER:: pi=3.14159265358979323846264338327950
	REAL(8), PARAMETER, DIMENSION(1:3)::weight=(/5./18.,4./9.,5./18./), &
	xq=(/-0.387298334620741688517926539978240, 0.0, &
	0.387298334620741688517926539978240/)
      CONTAINS
      !**************************************************************
      ! 5th order finite volume WENO reconstruction for two cell ends
      !**************************************************************
      subroutine weno_left(average,uleft)
      implicit none
      real, intent(in):: average(-2:2)
      real, intent(inout):: uleft
      real :: um(1:3), up(1:3), beta(1:3),alpha(1:3), ww(1:3), alpha_t(1:3), ww_t(1:3)
      integer:: i,k, is
      do k=1,3
         up(k)=0d0
         is=k-3
         do i=0,2
            up(k)=up(k)+c(-is-1,i)*average(is+i)
         enddo
      enddo
      
      beta(1)=13d0/12d0*(average(0)-2d0*average(1)+average(2))**2+0.25d0*(3d0*average(0)-4d0*average(1)+average(2))**2
      beta(2)=13d0/12d0*(average(-1)-2d0*average(0)+average(1))**2+0.25d0*(average(-1)-average(1))**2
      beta(3)=13d0/12d0*(average(-2)-2d0*average(-1)+average(0))**2+0.25d0*(average(-2)-4d0*average(-1)+3d0*average(0))**2
      
      alpha_t(1)=0.1/(1d-12+beta(1))**2
      alpha_t(2)=0.6/(1d-12+beta(2))**2
      alpha_t(3)=0.3/(1d-12+beta(3))**2
      
      ww_t(1)=alpha_t(1)/(alpha_t(1)+alpha_t(2)+alpha_t(3))   
      ww_t(2)=alpha_t(2)/(alpha_t(1)+alpha_t(2)+alpha_t(3))
      ww_t(3)=alpha_t(3)/(alpha_t(1)+alpha_t(2)+alpha_t(3))
      
      uleft=ww_t(1)*up(3)+ww_t(2)*up(2)+ww_t(3)*up(1)
      end subroutine weno_left
      
      
      
      !**************************************************************
      ! 5th order finite volume WENO reconstruction for two cell ends
      !**************************************************************
      subroutine weno_right(average,uright)
      implicit none
      real, intent(in):: average(-2:2)
      real, intent(inout):: uright
      real :: um(1:3), up(1:3), beta(1:3),alpha(1:3), ww(1:3), alpha_t(1:3), ww_t(1:3)
      integer:: i,k, is
      do k=1,3
         um(k)=0d0
         is=k-3
         do i=0,2
            um(k)=um(k)+c(-is,i)*average(is+i)
         enddo
      enddo
      
      beta(1)=13d0/12d0*(average(0)-2d0*average(1)+average(2))**2+0.25d0*(3d0*average(0)-4d0*average(1)+average(2))**2
      beta(2)=13d0/12d0*(average(-1)-2d0*average(0)+average(1))**2+0.25d0*(average(-1)-average(1))**2
      beta(3)=13d0/12d0*(average(-2)-2d0*average(-1)+average(0))**2+0.25d0*(average(-2)-4d0*average(-1)+3d0*average(0))**2
      
      alpha(1)=0.3/(1d-12+beta(1))**2
      alpha(2)=0.6/(1d-12+beta(2))**2
      alpha(3)=0.1/(1d-12+beta(3))**2
      ww(1)=alpha(1)/(alpha(1)+alpha(2)+alpha(3))   
      ww(2)=alpha(2)/(alpha(1)+alpha(2)+alpha(3))
      ww(3)=alpha(3)/(alpha(1)+alpha(2)+alpha(3))
      
      uright=ww(1)*um(3)+ww(2)*um(2)+ww(3)*um(1)
      end subroutine weno_right
            
      !**************************************************************
      ! 5th order finite volume WENO reconstruction for two cell ends
      !**************************************************************
      subroutine weno(average,uleft,uright)
      implicit none
      real, intent(in):: average(-2:2)
      real, intent(inout):: uleft, uright
      real :: um(1:3), up(1:3), beta(1:3),alpha(1:3), ww(1:3), alpha_t(1:3), ww_t(1:3)
      integer:: i,k, is
      do k=1,3
         um(k)=0d0
         up(k)=0d0
         is=k-3
         do i=0,2
            um(k)=um(k)+c(-is,i)*average(is+i)
            up(k)=up(k)+c(-is-1,i)*average(is+i)
         enddo
      enddo
      beta(1)=13d0/12d0*(average(0)-2d0*average(1)+average(2))**2+0.25d0*(3d0*average(0)-4d0*average(1)+average(2))**2
      beta(2)=13d0/12d0*(average(-1)-2d0*average(0)+average(1))**2+0.25d0*(average(-1)-average(1))**2
      beta(3)=13d0/12d0*(average(-2)-2d0*average(-1)+average(0))**2+0.25d0*(average(-2)-4d0*average(-1)+3d0*average(0))**2
      alpha(1)=0.3/(1d-12+beta(1))**2
      alpha(2)=0.6/(1d-12+beta(2))**2
      alpha(3)=0.1/(1d-12+beta(3))**2
      alpha_t(1)=0.1/(1d-12+beta(1))**2
      alpha_t(2)=0.6/(1d-12+beta(2))**2
      alpha_t(3)=0.3/(1d-12+beta(3))**2
      ww(1)=alpha(1)/(alpha(1)+alpha(2)+alpha(3))   
      ww(2)=alpha(2)/(alpha(1)+alpha(2)+alpha(3))
      ww(3)=alpha(3)/(alpha(1)+alpha(2)+alpha(3))
      ww_t(1)=alpha_t(1)/(alpha_t(1)+alpha_t(2)+alpha_t(3))   
      ww_t(2)=alpha_t(2)/(alpha_t(1)+alpha_t(2)+alpha_t(3))
      ww_t(3)=alpha_t(3)/(alpha_t(1)+alpha_t(2)+alpha_t(3))
      
      uright=ww(1)*um(3)+ww(2)*um(2)+ww(3)*um(1)
      uleft=ww_t(1)*up(3)+ww_t(2)*up(2)+ww_t(3)*up(1)
      end subroutine weno
      !***************************************************************************************
      ! 5th order finite volume WENO reconstruction for 3-point Gaussian quadrature points
      !***************************************************************************************
      subroutine weno_gauss(average,point_value)
      implicit none
      real, intent(in):: average(-2:2)
      real, intent(inout):: point_value(3)
      real :: um(1:3), up(1:3), beta(1:3),alpha(1:3), ww(1:3), alpha_t(1:3), ww_t(1:3)
      integer:: i,k, is, m
      beta(1)=13d0/12d0*(average(0)-2d0*average(1)+average(2))**2+0.25d0*(3d0*average(0)-4d0*average(1)+average(2))**2
      beta(2)=13d0/12d0*(average(-1)-2d0*average(0)+average(1))**2+0.25d0*(average(-1)-average(1))**2
      beta(3)=13d0/12d0*(average(-2)-2d0*average(-1)+average(0))**2+0.25d0*(average(-2)-4d0*average(-1)+3d0*average(0))**2
      DO m=1,3
         do k=1,3
            um(k)=0d0
            is=k-3
            do i=0,2
               um(k)=um(k)+cq(-is,i,m)*average(is+i)
            enddo
         enddo
	 if (m==2) then
	    !****splitting technique to remove the negative weight
	    alpha(1)=rq_p(1)/(1d-12+beta(1))**2
            alpha(2)=rq_p(2)/(1d-12+beta(2))**2
            alpha(3)=rq_p(3)/(1d-12+beta(3))**2
            ww(1)=alpha(1)/(alpha(1)+alpha(2)+alpha(3))   
            ww(2)=alpha(2)/(alpha(1)+alpha(2)+alpha(3))
            ww(3)=alpha(3)/(alpha(1)+alpha(2)+alpha(3))
	    alpha_t(1)=rq_n(1)/(1d-12+beta(1))**2
            alpha_t(2)=rq_n(2)/(1d-12+beta(2))**2
            alpha_t(3)=rq_n(3)/(1d-12+beta(3))**2
            ww_t(1)=alpha_t(1)/(alpha_t(1)+alpha_t(2)+alpha_t(3))   
            ww_t(2)=alpha_t(2)/(alpha_t(1)+alpha_t(2)+alpha_t(3))
            ww_t(3)=alpha_t(3)/(alpha_t(1)+alpha_t(2)+alpha_t(3))
	    point_value(m)=2.675*(ww(1)*um(3)+ww(2)*um(2)+ww(3)*um(1))-1.675*(ww_t(1)*um(3)+ww_t(2)*um(2)+ww_t(3)*um(1))
	 else
         alpha(1)=rq(1,m)/(1d-12+beta(1))**2
         alpha(2)=rq(2,m)/(1d-12+beta(2))**2
         alpha(3)=rq(3,m)/(1d-12+beta(3))**2
         ww(1)=alpha(1)/(alpha(1)+alpha(2)+alpha(3))   
         ww(2)=alpha(2)/(alpha(1)+alpha(2)+alpha(3))
         ww(3)=alpha(3)/(alpha(1)+alpha(2)+alpha(3))
         point_value(m)=ww(1)*um(3)+ww(2)*um(2)+ww(3)*um(1)
	 endif
      ENDDO
      end subroutine weno_gauss
      
      real function eval(aa,m,x0)
      implicit none 
      ! evaluate the value of the polynomial of degree m with coefficient a,
      ! at the location x=x0
      ! a(0)+a(1)*x0+...
      real, intent(in) :: aa(0:m), x0
      integer, intent(in) :: m
      integer :: i
    
      eval=aa(m)
      do 1 i=m-1,0,-1
           eval=eval*x0+aa(i)
     1 continue
      return
      end function eval
      
      END MODULE weno5