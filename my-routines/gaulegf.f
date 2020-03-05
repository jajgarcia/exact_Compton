! gauleg.f90     P145 Numerical Recipes in Fortran
! compute x(i) and w(i)  i=1,n  Legendre ordinates and weights
! on interval -1.0 to 1.0 (length is 2.0)
! use ordinates and weights for Gauss Legendre integration
!
      subroutine gaulegf(x1, x2, x, w, n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: x1, x2
      double precision, dimension(n), intent(out) :: x, w
      integer :: i, j, m
      double precision :: p1, p2, p3, pp, xl, xm, z, z1
      double precision, parameter :: eps=3.d-14
c      
      m = (n+1)/2
      xm = 0.5d0*(x2+x1)
      xl = 0.5d0*(x2-x1)
      do i=1,m
         z = cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
         z1 = 0.0
        do while(abs(z-z1) .gt. eps)
           p1 = 1.0d0
           p2 = 0.0d0
           do j=1,n
              p3 = p2
              p2 = p1
              p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
           end do
           pp = n*(z*p1-p2)/(z*z-1.0d0)
           z1 = z
           z = z1 - p1/pp
        end do
        x(i) = xm - xl*z
        x(n+1-i) = xm + xl*z
        w(i) = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
        w(n+1-i) = w(i)
      end do
      end subroutine gaulegf
