      program drive_SRF
      use omp_lib
c     Driving program for the writting of the Super Redistribution Function SRF
c     due to Compton scattering (see more details in super_Compton_RF.f).
c     The output of this code is used by the XILLVER code.
c
c     It writes to an ascii file for the moment.
c
c     To-Do:
c    - The minimum number of angles (mgi) needed for accurate angular integration
c      depends on temp and initial photon energy. 3000 seems to be the worse case
c      scenario, but we could make mgi to change according to the conditions for speed up.
c
c     Authors: Javier Garcia (javier@caltech.edu)
c              Ekaterina Sokolova-Lapa (ekaterina.sokolova-lapa@fau.de)
c
c........
      implicit none
      integer nmaxp, itrans, mgi, ii

      parameter (nmaxp=5000, itrans=70, mgi=3000)
      double precision pemin, pemax, pemax2
      double precision theta(itrans), wp(nmaxp), df(nmaxp)
      double precision skn(nmaxp,itrans)
      double precision ikbol, ergsev, mec2, pi
      double precision smit(mgi), agt(mgi)
      double precision tini, tfin, tcpu, temp
c
c     Get current time
      call cpu_time(tini)
c
c
      ergsev  = 1.602197d-12     ! Convert eV to ergs
      ikbol   = 1.16d4           ! inverse of kbol (K * ev-1)
      pi      = 4.d0*datan(1.d0) ! pi number
      mec2    = 5.11d5           ! m_e c^2 (eV)
c
c     Array of temperatures
      temp = 1.d4                ! Gas temp in K
      do ii=1,itrans
         theta(ii) = temp/ikbol/mec2
         temp = temp*(1.d10/1.d4)**(1.d0/dfloat(itrans-1))
      enddo
c
c     Photon energy grid
      pemin = 0.1d0
      pemax = 9.9d5
      pemax2 = 1.d6
      call enegrd(nmaxp, pemin, pemax, pemax2,
     1            wp, df)
c
c     Calculate the Compton Cross Section
      call scattxs(nmaxp, wp, itrans, theta,
     1             skn)
c
c     Get the Gaussian quadratures for angular integration
      call gaulegf(-1.d0, 1.d0, smit, agt, mgi)
c      print *,smit
c      smit1(1)=-1.d0
c      smit1(mgi+2)=1.d0
c      agt1(1)=smit(1)-(-1.d0)
c      agt1(mgi+2)=1.d0-smit(mgi)
c      do ii=1,mgi
c            smit1(ii+1)=smit(ii)
c            agt1(ii+1)=agt(ii)
c      enddo
c      agt1(2)=agt1(2)-agt1(1)
c      agt1(mgi+1)=agt1(mgi+1)- agt1(mgi+2)
c      print*,agt1
c      print*,smit1
c     Produce file with all SRF's
      call super_Compton_RF_fits(itrans, theta, nmaxp, wp, df, skn,
     1     mgi, smit, agt)
c
c     Get current time
      call cpu_time(tfin)
      tcpu=tfin-tini
c
      print *, ' '
      print *, 'CPU time (s) =',real(tcpu)
c
      end program drive_SRF
