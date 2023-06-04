.. _tools_panglos:

=======
Panglos
=======

Introduction
------------

Panglos is a processing code for thermal scattering laws. It reads an ENDF-6
formatted file with a File 7, and generates a custom ACE file which can be read
by PapillonNDL. It has several advantages over the more traditional THERMR
processing in NJOY.

1. For Incoherent Inelastic scattering, the cross section is initially
   calculated on the traditional THERMR energy grid, but is then linearized,
   adding points where necessary. This allows finer details in the cross section
   to be exposed, such as for H1 in ZrH.

2. For Incoherent Inelastic scattering, the Direct S(a,b) method, propossed by
   Ballinger is used to model the scattering distribution. This allows for a
   continuous distribution in both energy and angle. Currently, THERMR only
   provides continuous energy distributions, and uses discrete angles.

3. For Incoherent Elastic scattering, only the characteristic bound cross
   section and the Debye Waller integral are provided. This allows use of the
   analytic equations for both the integral cross section, and the angular
   distribution.

IncoherentInelastic
-------------------

.. doxygenclass:: IncoherentInelastic

IncoherentElastic
-----------------

.. doxygenclass:: IncoherentElastic

CoherentElastic
---------------

.. doxygenclass:: CoherentElastic

Sab
---

.. doxygenclass:: Sab

ShortCollisionTimeSab
---------------------

.. doxygenclass:: ShortCollisionTimeSab

FreeGasSab
---------------------

.. doxygenclass:: FreeGasSab

Tabulatedsab
------------

.. doxygenclass:: TabulatedSab

ENDFSab
-------

.. doxygenclass:: ENDFSab