!-------------------------------------------------------------------------------------------------
      subroutine super_Compton_RF_fits(itrans, theta, nmaxp, wp, df,
     & skn,  mgi, smit, agt)
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
      real*8 check1(nmaxp), limit, pmin, pmax(nmaxp)
      real*8 prob(nmaxp,nmaxp), srf
      real*8 mec2, ecen, isp, ikbol, temp(itrans),check
      integer jj, kk, point(nmaxp), indices(nmaxp,nmaxp),ind_arr(nmaxp)
      integer :: n ,indmax(nmaxp)
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


      do iz=1, itrans
            temp(iz) = theta(iz)*ikbol*mec2
c            print *,temp(iz)
      enddo
!crate and fill the first extension with temperature and energy
      call write_param(itrans, nmaxp, temp, wp, filename)

!append and other extension (3rd) to store the SRF
!IMPORTANT: this routine leaves the fits file opened      
      call add_HDU(itrans, nmaxp, filename, unit)
      do iz = 1, itrans
         do kk=1,nmaxp
            point(kk)=0
            pmax(kk)=0.d0
            check1(kk) = 0.d0
            indmax(kk)=0
            do np = 1, nmaxp
               prob(np,kk) = 0.d0
               indices(np,kk)=0
            enddo
         enddo
c
c         temp = theta(iz)*ikbol*mec2          ! temperature in K
         do np = 1, nmaxp
            ecen = wp(np)
c$omp parallel num_threads(30)
c$omp& shared(iz,wp,prob,point,pmax,nmaxp,np,ecen,temp,mgi,smit,
c$omp& agt,mec2,indmax)
c$omp& private(jj)
c$omp do
            do jj=1,nmaxp
               call probab(temp(iz),wp(jj)/mec2,ecen/mec2,mgi,smit
     1          ,agt,prob(jj,np))
               if(prob(jj,np).gt.pmax(np)) then
                  pmax(np) = prob(jj,np)
                  indmax(np)=jj
               endif
            enddo
c$omp end do
c$omp end parallel
         if(pmax(np).le.0.d0) then
            print *,np
         endif
         enddo
         do np=1,nmaxp
            check=0.d0
            do jj=indmax(np),1,-1
                  kk=point(jj)+1
                  if(prob(jj,np).ge.(pmax(np)*limit))then
                        indices(jj,kk)=np
                        check=check+df(jj)*prob(jj,np)
                        point(jj)=kk
                  else
                        exit
                  endif
            enddo

            do jj=indmax(np)+1,nmaxp
                  kk=point(jj)+1
                  if(prob(jj,np).ge.(pmax(np)*limit))then
                        indices(jj,kk)=np
                        check=check+df(jj)*prob(jj,np)
                        point(jj)=kk
                  else
                        exit
                  endif

            enddo
            check1(np)=check
         enddo

         
!     From here crated by Gullo Dec 2020
!     Once it calculates the srf it writes directly the fits file
!     It needs to differentiate the first call, where it writes the extension from all the other calls 
      do jj = 1, nmaxp
c         if (jj.lt.300)then
c         print *,"*****",jj
c         do kk=1,point(jj)
c            print*,indices(jj,kk)
c         enddo
c         endif
         do kk=1,nmaxp
            srf_arr(kk) = 0.d0
c            ind_arr(kk) = 0
         enddo
         do kk = 1, point(jj)
            np=indices(jj,kk)
            srf = prob(jj,np)*skn(np,iz)*df(np)/wp(np)/check1(np)
!save the SRF array in order to pass it to the fits file routine   
            srf_arr(kk) = srf
c            ind_arr(kk) = np
         enddo
         call add_row_HDU(n, nmaxp, point(jj), indices(jj,1),
     &        skn(jj,iz), srf_arr, unit)
            n = n + 1           ! increment the row number
      enddo
               
      print *,iz
      enddo      
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
