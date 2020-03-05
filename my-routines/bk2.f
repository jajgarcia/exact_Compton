c........
      double precision function bk2 (x)
      implicit none
      double precision i0,i1,k0,k1,x,t,z
c.......  Computations of the modified Bessel function of the second
c.......  order, scaled by exponent
c.......
c.......             K_2 (x) \times exp(x)  -  for  X.GT.0,
c.......
c.......  Method: asymptotic expansions, given by Abramowitz and Stegun
c.......  are used for computing I_0, I_1, K_0, and K_1.
c.......
c.......  Final step makes use of recurrence relation:
c.......         K_{n+1} (x) = {2n/x} K_n (x) + K_{n-1} (x)
c.......
c.......  The function has been tested in the range X: between 0.01 and 500.
c.......  (Madej 1989).
c.......
c
c.....High temperature (kT > 255.5 keV)
c      if (2.d0-x.ge.0.d0)then
      if (x.le.2.d0) then
         z=(x/2.d0)**2
         t=(x/3.75d0)**2
         i0 = (((((0.0045813d0*t+0.0360768d0)*t+0.2659732d0)*t
     *        +1.2067492d0)*t+3.0899424d0)*t+3.5156229d0)*t + 1.d0
         i1 = (((((.00032411d0*t+.00301532d0)*t+.02658733d0)*t
     *        +.15084934d0)*t+.51498869d0)*t+.87890594d0)*t + 0.5d0
         i1=x*i1
         k0 = (((((.00000740d0*z+.00010750d0)*z+.00262698d0)*z
     *         +.03488590d0)*z+.23069756d0)*z+.42278420d0)*z
     *         -.57721566d0 - i0*dlog(x/2.d0)
         k1 = (((((-.00004686d0*z-.00110404d0)*z-.01919402d0)*z
     *         -.18156897d0)*z-.67278579d0)*z+.15443144d0)*z  + 1.d0
         k1 = k1/x + i1*dlog(x/2.d0)
c.......         K_{n+1} (x) = {2n/x} K_n (x) + K_{n-1} (x)
         z = k0 + 2.d0/x*k1
         bk2 = z * dexp(x)
c
c.....This is within the limit this Bessel function was tested
      else if (x.lt.500.d0) then
         z=2.d0/x
         k0 = (((((0.00053208d0*z-0.00251540d0)*z+0.00587872d0)*z
     *        -0.01062446d0)*z+0.02189568d0)*z-0.07832358d0)*z
     *        +1.25331414d0
         k1 = (((((-0.00068245d0*z+0.00325614d0)*z-0.00780353d0)*z
     *        +0.01504268d0)*z-0.03655620d0)*z+0.23498619d0)*z
     *        +1.25331414d0
         bk2 = (k0+z*k1)/dsqrt(x)
c
c.....For very low temp (< 10^7 K). Proposed by Katya, but so far seems to
c.....provide the same answer (so maybe we don't need it).
      else
         z=1.d0/x
         bk2 = dsqrt(3.14159265d0/2.d0)*(dsqrt(z)+1.875d0*z**(1.5d0)
     *         +0.8203125d0*z**(2.5d0))
      endif
      return
      end
c

