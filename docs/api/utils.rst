.. _api_misc:

=========
Utilities
=========

NDLibrary
-----------

.. doxygenclass:: pndl::NDLibrary


MCNPLibrary
-----------

.. doxygenclass:: pndl::MCNPLibrary

SerpentLibrary
--------------

.. doxygenclass:: pndl::SerpentLibrary

ACE
---

.. doxygenclass:: pndl::ACE

XSPacket
--------

.. doxygenstruct:: pndl::XSPacket

PCTable
-------

.. doxygenclass:: pndl::PCTable

Frame
-----

.. doxygenenum:: pndl::Frame

CMToLab
-------

.. doxygenstruct:: pndl::CMToLab

Interpolation
-------------

.. doxygenenum:: pndl::Interpolation

Interpolator
------------

.. doxygenclass:: pndl::Interpolator

Linearization
-------------

.. doxygenfunction:: pndl::linearize(const std::vector<double> &x, const std::vector<double> &y, std::function<double(double)> f, double tolerance = 0.001)

.. doxygenfunction:: pndl::linearize(double x_min, double x_max, std::function<double(double)> f, double tolerance = 0.001)

Random Number Generation
------------------------

.. doxygenfunction:: pndl::rng

.. doxygenfunction:: pndl::rng_seed

.. doxygenfunction:: pndl::rng_reset

.. doxygenfunction:: pndl::rng_advance

ZAID
----

.. doxygenclass:: pndl::ZAID

Element
-------

.. doxygenclass:: pndl::Element

Isotope
-------

.. doxygenclass:: pndl::Isotope

Nuclide
-------

.. doxygenclass:: pndl::Nuclide

PNDLException
-------------

.. doxygenclass:: pndl::PNDLException
