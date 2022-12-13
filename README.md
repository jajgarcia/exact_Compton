# exact_Compton
This is version 0.3.0 of the exact_Compton model.

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
