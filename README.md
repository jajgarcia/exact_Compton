# exact_Compton
These routines compute the exact redistribution function due to Compton scattering

Included is a driving program (driveSRF.f) which calls all pertinent routines.

The main output of the code is what we calle the "super redistribution function
(SRF)" for Compton scatterting. For a given gas temperature T and final photon
energy Ef, the SRF is defined for a set of initial photon energies Ei as:

         SRF(T,Ef,Ei) = IRF(Ef,Ei)/N(Ei)*skn(Ei)*dEi/Ei
     
where IRF(Ef,Ei) is in fact the inverse redistribution function for the Compton
scattering of a photon from initial energy Ei to final energy Ef; N(Ei) is the
normalization to ensure photon number conservation; and skn(Ei) is the
Klein-Nishina cross section.  This routine implements the exact Compton RF from
Madej et al. (2017).  The SRF contains all the information needed for the
convolution of a given spectrum to account for the Compton scattering at the
given temperature.  Only significant values of the SRF are actually written,
i.e., when RF > limit.

The output of this code is used by the XILLVER model.

     Version: 0.4.0 - Wed Nov 23rd 14:15:39 PDT 2022

     Authors: Javier Garcia (javier@caltech.edu)
              Ekaterina Sokolova-Lapa (ekaterina.sokolova-lapa@fau.de)
                (see Garcia et al. 2020 in prep)
	      Jameson Dong (jdong2@caltech.edu)
	      Isabel Franco (francog@caltech.edu)
	      Guglielmo Mastroserio (guglielmo.mastroserio@inaf.it)
	      With routines provided by J. Madej and A. Rozanska (see Madej et al. 2017).
