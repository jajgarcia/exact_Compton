      subroutine crsexact(theta, eps, totalcrs)
c........
c........ Total Compton scattering cross section in a gas of
c........ relativistic electrons. 
c........ Given by Poutanen & Svensson, 1996

c........ electron temperature ('temper') and energy ('en')
c........ are in mc^2
c........ totalcrs is returned in [cm^-2]

c........ NOTE:
c........ test in range of electron temperature [8e-3 - 1.7] mc^2
c........ (around [5e+7 - 1e+10] Kelvins).
c........ below 1e+8 K the Klein-Nishina cross section can be used.


      implicit none
      integer i, gi, ki
      parameter (gi = 200)
      double precision eps, theta, totalcrs
      double precision mec2, sigma_t, g_m1, aj, mn, x, bk2, crsgm
      parameter (sigma_t = 6.65246d-25)
      dimension g_m1(gi), aj(gi)

      if (eps.le.1d-4) then
         totalcrs = sigma_t
      else
         totalcrs = 0d0

         mn = 3.d0 * sigma_t / (16.d0 * eps * eps * theta)

c     integration over gamma

         call gaulegf(1d0, 20d0, g_m1, aj, gi)
         do i = 1,gi
            totalcrs = totalcrs + aj(i)*crsgm(theta, eps, g_m1(i))
         enddo

         x = 1d0/theta
         totalcrs = mn * totalcrs / bk2(x)

      endif
      !write (2,'(1p2e15.5)') eps, totalcrs

      return
      end
      double precision function crsgm(theta, eps, gm)
c     integrand of the total exact relativistic cross section
c     for Compton scattering (for each gamma point)      
      implicit none
      integer i, nksi
      parameter ( nksi = 20 )
      double precision ksi(nksi), wksi(nksi), x
      double precision theta, eps, integ, gm, fep, z, a, b, bk2

      crsgm = 0d0

      z = dsqrt(gm*gm - 1d0)

c     boarders and point for internal integral over ksi: 
      a = eps*(gm - z)
      b = eps*(gm + z)
      call gaulegf(a, b, ksi, wksi, nksi)
c     integral over ksi: 
      integ = 0d0
      do i =1,nksi
         integ = integ + wksi(i)*dlog(1d0 + 2d0*ksi(i))/ksi(i)
      enddo

      crsgm = ( (eps*gm + 4.5d0 + 2d0*gm/eps)
     **dlog((1d0 + 2d0*eps*(gm+z)) / (1d0 + 2d0*eps*(gm-z)))
     *-2d0*eps*z + z*(eps-2d0/eps)
     **dlog(1d0+4d0*eps*gm + 4d0*eps*eps)
     *+ (4d0*eps*eps*z*(gm + eps))
     */(1d0 + 4d0*eps*gm + 4d0*eps*eps) - 2d0*integ)

c     exponent from relativistic Maxwellian distribution
      x = 1d0/theta
      fep = dexp((1.-gm)*x)

      crsgm = crsgm * fep

      return
      end


