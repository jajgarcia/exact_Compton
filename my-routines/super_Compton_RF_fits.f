!-------------------------------------------------------------------------------------------------
      subroutine super_Compton_RF_fits(itrans, theta, nmaxp, wp, df,
     1                                 skn,  mgi, smit, agt)
c
c     This routine writes a file with the super redistribution function (SRF) for Compton
c     scatterting. For a given gas temperature T and final photon energy Ef, the SRF is
c     defined for a set of initial photon energies Ei as:
c
c         SRF(T,Ef,Ei) = IRF(Ef,Ei)/N(Ei)*skn(Ei)*dEi/Ei
c
c     where IRF(Ef,Ei) is in fact the inverse redistribution function for the Compton scattering
c     of a photon from initial energy Ei to final energy Ef; N(Ei) is the normalization to
c     ensure photon number conservation; and skn(Ei) is the Klein-Nishina cross section.
c     This routine implements the exact Compton RF from Madej et al. (2017).
c     The SRF contains all the information needed for the convolution of a given spectrum
c     to account for the Compton scattering at the given temperature.
c     Only significant values of the SRF are actually written, i.e., when RF > limit.
c
c     Input arguments:
c         itrans: Total number of temperatures
c         theta: Array (itrans) of temperatures in kT/mec2
c         nmaxp: Total number of energy points
c         wp: Array (nmaxp) of energies in eV
c         df: Array (nmaxp) of delta energies in eV
c         skn: Array (nmaxp) of Kein-Nishina cross sections
c         mgi: Total number of angles
c         smit: Array (mgi) Legendre ordinates (angles)
c         agt: Array (mgi) weights
c
c     Output arguments:
c         None
c
c     Requires:
c         probab.f: Routine for the RF calculation
c
      implicit none
      integer itrans, nmaxp, mgi, iz, np
      real*8 theta(itrans), wp(nmaxp), df(nmaxp), skn(nmaxp,itrans)
      real*8 smit(mgi), agt(mgi)
      real*8 check1(nmaxp), limit, pmin
      real*8 prob(nmaxp,nmaxp), srf
      real*8 mec2, ecen, isp, ikbol, temp
      integer jj, kk, point(nmaxp), indices(nmaxp,nmaxp)
      integer :: n
      character (len=200) filename
! Added by Gullo
c$$$  double precision, allocatable :: srf_arr(:)
      integer          :: unit, status
      double precision :: srf_arr(nmaxp)
c Initialize status
      status=0
c
      ikbol   = 1.16d4      ! inverse of kbol (K * ev-1)
      mec2  = 5.11d5        ! m_e c^2 (eV)
      isp = 0.5641895835d0  ! 1/sqrt(pi)
      limit = 1.d-3         ! Limit for the redistribution function
c
c$$$100   format(2i8)
c$$$101   format(i8,ES16.8E3)
c$$$102   format(i8,ES16.8E3,i8,ES16.8E3)
c$$$103   format(2i8,2ES16.8E3)
c
c     Check1 is to ensure photon number is conserved in scatterings

      filename = 'table.fits' !name of the fits file
      n = 1 ! column number of the fits file

      call create_fits(filename) !create the fits file

!crate and fill the first extension with temperature and energy
      call write_param(itrans, nmaxp, theta, wp, filename)

!append and other extension (3rd) to store the SRF
!IMPORTANT: this routine leaves the fits file opened
      call add_HDU(itrans, nmaxp, filename, unit)

      do 1001 iz = 1, itrans
         do kk=1,nmaxp
            point(kk)=0
            do np = 1, nmaxp
               prob(np,kk) = 0.d0
            enddo
         enddo
c
         temp = theta(iz)*ikbol*mec2          ! temperature in K
c
         do 989 np = 1, nmaxp
c
            check1(np) = 0.d0
c           Approx energy of the probability maximum
            ecen = wp(np)
c           Probability at ecen
            call probab(temp,ecen/mec2,ecen/mec2,mgi,smit,agt,
     1                  prob(np,np))
            check1(np) = check1(np) + df(np) * prob(np,np)              ! Normalization
            pmin = prob(np,np)*limit
            kk = point(np)+1
            indices(np,kk)=np
            point(np)=kk
c
c           Integration to the left
            do jj = np-1, 1, -1
               kk=point(jj)+1
               call probab(temp,wp(jj)/mec2,wp(np)/mec2,mgi,smit,agt,
     1                     prob(jj,np))
               if (prob(jj,np).ge.pmin)then                             ! Stop at pmin
                   check1(np) = check1(np) + df(jj) * prob(jj,np)       ! Normalization
                   indices(jj,kk)=np
                   point(jj)=kk
               else
                   exit
               endif
            enddo
c
c           Integration to the right
            do jj = np+1, nmaxp
               kk=point(jj)+1
               call probab(temp,wp(jj)/mec2,wp(np)/mec2,mgi,smit,agt,
     1                     prob(jj,np))
               if (prob(jj,np).ge.pmin)then                             ! Stop at pmin
                   check1(np) = check1(np) + df(jj) * prob(jj,np)       ! Normalization
                   Indices(jj,kk)=np
                   point(jj)=kk
               else
                   exit
               endif
            enddo
989      enddo


!     From here crated by Gullo Dec 2020
!     Once it calculates the srf it writes directly the fits file
!     It needs to differentiate the first call, where it writes the extension from all the other calls
      do jj = 1, nmaxp
         do kk = 1, point(jj)
            np=indices(jj,kk)
            srf = prob(jj,np)*skn(np,iz)*df(np)/wp(np)/check1(np)
!save the SRF array in order to pass it to the fits file routine
            srf_arr(kk) = srf
         enddo
            call add_row_HDU(n, nmaxp, point(jj), indices(jj,1),
     &        skn(jj,iz), srf_arr, unit)
            n = n + 1           ! increment the row number
      enddo

c
1001  continue
c

c The FITS file must always be closed before exiting the program.
c Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      write(*,*) 'final status', status
      call ftclos(unit, status)
      write(*,*) 'final status', status
      call ftfiou(unit, status)
      write(*,*) 'final status', status
      if (status .gt. 0) then
         call printerror(status)
c$$$         write(*,*) 'end file'
      endif
      return
      end subroutine
