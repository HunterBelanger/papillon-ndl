# Panglos

Panglos is a small code which processes thermal scattering laws from MF 7 in the
ENDF-6 format. It produces a customized ACE file, which contains all of the
formation required to sample continuous distributions for thermal scattering
laws. This includes:

- The Debye-Waller integral for sampling Incoherent Elastic scattering in an
  exact manner.
  
- The continuous distributions to sample the S(a,b) for the energy and
  scattering angle continuously, for Incoherent Inelastic scattering. This data
  processing is based on the direct S(a,b) method of Ballinger.
