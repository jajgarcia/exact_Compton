      subroutine enegrd(nmaxp, ppemin, ppemax, ppemax2,  !inp
     1                  wp, df)                          !out
      implicit none
      integer nmaxp
      real*8 wp(nmaxp), ppemin, ppemax, ppemax2
      integer numcon, numcon2, numcon3, ll, ll2
      real*8 ebnd1, ebnd2, ebnd2o, dele, df(nmaxp)
c
      numcon = nmaxp
      if (numcon.lt.4) stop 'in ener: numcon error'
      numcon2=max(2,nmaxp/50)
      numcon3=numcon-numcon2
      ebnd1=ppemin
      ebnd2=ppemax
      ebnd2o=ebnd2
      dele=(ebnd2/ebnd1)**(1./dfloat(numcon3-1))
      wp(1)=ebnd1
      do ll=2,numcon3
        wp(ll)=wp(ll-1)*dele
      enddo
      ebnd2=ppemax2
      ebnd1=ebnd2o
      dele=(ebnd2/ebnd1)**(1./dfloat(numcon2-1))
      do ll2=1,numcon2
        ll=ll2+numcon3
        wp(ll)=wp(ll-1)*dele
      enddo
c
      df(1)=0.d0
      do ll=2,nmaxp
         df(ll) = wp(ll) - wp(ll-1)
      enddo
      df(1)=df(2)  !! Need to check this!
c
      return
      end
