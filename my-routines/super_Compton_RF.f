      subroutine super_Compton_RF(itrans, theta, nmaxp, wp, df, skn,  !inp
     1                            mgi, smit, agt)                     !inp
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
c
      ikbol   = 1.16d4      ! inverse of kbol (K * ev-1)
      mec2  = 5.11d5        ! m_e c^2 (eV)
      isp = 0.5641895835d0  ! 1/sqrt(pi)
      limit = 1.d-3         ! Limit for the redistribution function
c
100   format(2i8)
101   format(i8,ES16.8E3)
102   format(i8,ES16.8E3,i8,ES16.8E3)
103   format(2i8,2ES16.8E3)
c
c     Write the first line in output: number of temps and final energies
      open(10, status = 'unknown', file = 'srf.txt')               ! Set of SRF's for all temps
      write(10,100)itrans, nmaxp
c
c     Check1 is to ensure photon number is conserved in scatterings
c
      do 1001 iz = 1, itrans
         do kk=1,nmaxp
            point(kk)=0
            do np = 1, nmaxp
               prob(np,kk) = 0.d0
            enddo
         enddo
c
         temp = theta(iz)*ikbol*mec2                              ! temperature in K
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
                   indices(jj,kk)=np
                   point(jj)=kk
               else
                   exit
               endif
            enddo
989      enddo
c
c
c
c        Write down the SRF's for each temperature and final energy
c        On an ascii file at the moment, it should be into a FITS file
c        Also, we don't really need all the indices, but it's ok for now
c        Now also writes the precalculated scattering cross section skn
c
         write(10,101)iz,temp
c
         do jj = 1, nmaxp
             write(10,102)jj,wp(jj),point(jj),skn(jj,iz)
             do kk = 1, point(jj)
               np=indices(jj,kk)
               srf = prob(jj,np)*skn(np,iz)*df(np)/wp(np)/check1(np)
               write(10,103)kk,np,wp(np),srf
             enddo
         enddo
c
1001  continue
c
      close(10)
      return
      end subroutine
