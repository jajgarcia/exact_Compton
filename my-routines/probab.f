      subroutine probab(temper,enprim,encent,mgi,smit,agt,total)
c........ 
c........ Compton scattering in gas of relativistic thermal electrons.
c........
c........ Subroutine for the computations of:
c........   - angle-averaged Compton scattering probabilities
c........
c........ Algorithm used in computations:
c........  nr=1 : exact by Suleimanov et al. (2012), Madej et al. (2017)
c........
      implicit none
      integer i,idiv,mgi,in,j
      double precision smit(mgi),agt(mgi),xprime,step,xlow,xhigh
      double precision x,xdim,z,eps,eps1,temper,steplg
      double precision profil,redist1,total1
      double precision hhh,skkk,ccc,smhy,smel,esu,eee,sigma,pi,eve,evf
      double precision total,disp,ecentr,gauss,enprim,encent
      parameter (hhh= 6.62620d-27, skkk= 1.38062d-16, ccc= 2.99793d+10,
     *          smhy= 1.67333d-24, smel= 9.10956d-28, esu= 4.80325d-10,
     *           eee= 4.80325d-10, sigma=5.66961d-05, pi=3.141592654d0)
c........	       eV - associated energy	  erg
      parameter ( eve=1.60219d-12,  evf=2.41838d+14 )
      x=smel*ccc**2/skkk/temper
      eps=encent
      eps1=enprim
      xdim=eps*smel*ccc**2/hhh              ! in Hz      frequency
      xprime=eps1*smel*ccc**2/hhh           ! in Hz      frequency
c.......
c....... Computation of the redistribution probability (\nu,\nu^\prime,x)
c....... for the initial frequency \nu. Inverted temperature of electrons
c.......               x = m_e c**2 / kT 

c....... Integration over outgoing cosines:
      redist1=0.d0
      do j=1,mgi
         redist1 = redist1 + agt(j)*profil(1,eps,eps1,smit(j),x)
      enddo
      redist1 = 2.d0*pi*redist1     !*eps/xdim
c........
c........  This is the exact redistribution function in units sec,
c........  as in Madej et al. (2017) or Sazonov & Sunyaev (2000).
c........
c      write (3,'(1hc/1hc,6x,26hCOMPTON SCATTERING PROFILE/
c      .   1hc,4x,6hTEMP =,1pe10.3,13h K and FREQ =,1pe12.5,3h Hz,15x,
c      .   5hmgi =,i7/1hc)') temper,xdim,mgi
c      write (3,'(1hc,15x,10hFREQ PRIME,4x,11hPROBABILITY)')
c      write (3,'(1hc,17x,5hin Hz,10x,6hfor Hz/1hc)')
c      write (3,'(1hc,16x,7hE prime,11x,1hE,12x,5hexact,10x,6hJavier)')
c
c      in=1
c      disp=eps*dsqrt(2*skkk*temper/(smel*ccc**2)+0.4*eps**2)
c      ecentr=eps*(1+4*skkk*temper/(smel*ccc**2)-eps)
c      gauss=1.d0/(dsqrt(pi)*disp)*dexp(-(eps1-eps)**2/disp**2)
c      write (3,'(1p3e15.5)') eps1,eps,redist1
      total=redist1
c.......  The above gaussian was computed as in Garcia et al. (2013)
c.......  and renormalized to sec. It is not fully correct as above, 
c..,....  since still it must be divided by the Klein-Nishina opacity,
c.......  and now I do not remember their equation. 
      return
      end
c.......
c.......  nr=1 : exact by Suleimanov et al. (2012), Madej et al. (2017)
c.......
      double precision function profil(nr,eps,eps1,costh,x)
      implicit none
c.......  Computation of the exact Compton redisribution function
c.......	     P (epsilon, epsilon1, costh, x),
c.......  for photons scattered by free electrons in thermal gas.
c.......
c.......  epsilon = initial energy (in units of mc**2),
c.......  epsilon_1 = final energy of a scattered photon,
c.......  costh    = cosine of scattering angle, and
c.......
c.......	      x    = m c**2 / kT
c.......
c.......  Thirty-two point Gauss-Laguerre quadrature is used below.
c.......
      integer i,j,nr,nsize
      double precision integ,eps,eps1,costh,x,ag,aj,xj,bk2,gammin
      double precision pi,zfe,gam_star,q_up,q_down,fexact
      parameter ( nsize = 32 )
      dimension xj(nsize),aj(nsize)
      data pi/3.141592654d0/
c.......  The following gives nodes and weights of the 32-point
c.......  Gauss-Laguerre quadrature for integration over electron
c.......  Lorenz factor.
c.......  tables: XJ - nodes, AJ - weights.
      data (xj(i),aj(i),i=1,nsize)/0.044489366d0,1.09218342d-1,0.2345261
     *1d0,2.10443108d-1,0.576884629d0,2.352132297d-1,1.07244875d0,1.9590
     *3336d-1,1.72240878d0,1.29983786d-1,2.528336706d0,7.05786239d-2,3.4
     *92213273d0,3.17609125d-2,4.61645677d0,1.191821484d-2,5.9039585d0,3
     *.738816295d-3,7.358126733d0,9.8080331d-4,8.98294092d0,2.14864919d-
     *4,10.78301863d0,3.920342d-5,12.76369799d0,5.9345416d-6,14.93113976
     *d0,7.4164046d-7,17.29245434d0,7.6045679d-8,19.85586094d0,6.3506022
     *d-9,22.63088901d0,4.28138297d-10,25.62863602d0,2.30589949d-11,28.8
     *6210182d0,9.7993793d-13,32.34662915d0,3.2378017d-14,36.10049481d0,
     *8.17182344d-16,40.1457198d0,1.54213383d-17,44.509208d0,2.11979229d
     *-19,49.22439499d0,2.05442967d-21,54.3337213d0,1.34698259d-23,59.89
     *250916d0,5.6612941d-26,65.97537729d0,1.41856055d-28,72.68762809d0,
     *1.91337549d-31,80.18744698d0,1.19224876d-34,88.73534042d0,2.671511
     *22d-38,98.82954287d0,1.33861694d-42,111.7513981d0,4.5105362d-48/
c.......
      profil=0.d0
      q_down=eps*eps1*(1-costh)
      q_up=dsqrt((eps-eps1)**2+2*q_down)
      gam_star=(eps-eps1+q_up*dsqrt(1+2.d0/q_down))/2.d0
      integ=0.d0
c.......
c.......  Integration over gamma (Lorenz factor)
   11 do 16 j=1,nsize
         zfe=fexact(eps,eps1,costh,xj(j)/x+gam_star)
         integ=integ+zfe*aj(j)
c.......
c.......  In the following the final expression corresponds to Eq.13 of
c.......  Madej et al. (2017), where the profile is multiplied by the
c.......  factor eps1/eps*exp(eps-eps1) related to the detailed balance
c.......  condition.
c.......  Note: the unity in the exponent below accounts for the Bessel
c.......  function being scaled by exp(x). See bk2.f for the details.
   16 continue
      integ=integ*dexp((1-gam_star+eps-eps1)*x)
      profil=3.d0/32.d0/pi/bk2(x)*integ     *eps1/eps
c.......
      return
      end
c.......
c.......  Formulae from Madej et al. (2017):
      double precision function fexact (eps,eps1,costh,gamma)
      implicit none
      double precision eps,eps1,costh,gamma,a2_plus,a2_minus,
     *   q_up,q_down,d_plus,d_minus,a_minus,a_plus,a_diff,a2_diff
c.......
      q_down=eps*eps1*(1-costh)
      q_up=dsqrt((eps-eps1)**2+2*q_down)
      a2_minus=(gamma-eps)**2+(1+costh)/(1-costh)
      a2_plus=(gamma+eps1)**2+(1+costh)/(1-costh)
      d_minus=(a2_plus-a2_minus-q_up**2)/2.d0
      d_plus =(a2_plus-a2_minus+q_up**2)/2.d0
c.......
      a_minus=dsqrt(a2_minus)
      a_plus=dsqrt(a2_plus)
      a_diff=(-2*gamma*(eps+eps1)+eps**2-eps1**2)/(a_minus+a_plus)
      a2_diff=-2*gamma*(eps+eps1)+eps**2-eps1**2
c........ Final formula from Madej et al. (2017):
      fexact= 2/q_up-(q_down-2)/(q_down*a_plus*a_minus)*a_diff
     *  +((a2_minus+a_minus*a_plus+a2_plus)*q_up**2-a2_diff**2
     *  -a_minus*a_plus*a_diff**2)/(2*q_down**2
     *  *a_plus**3*a_minus**3)*a_diff
c........ Formula from Suleimanov et al. (2012):
c      fexact= 2/q_up + (q_down**2-2*q_down-2)/q_down**2
c     *     *(1/dsqrt(a2_minus)-1/dsqrt(a2_plus))
c     *     +(d_minus/a2_minus**1.5d0+d_plus/a2_plus**1.5d0)/q_down**2
      return
      end
c........
