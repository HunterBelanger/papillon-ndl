.. Papillon Nuclear Data Library documentation master file, created by
   sphinx-quickstart on Fri Feb 19 20:55:25 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Papillon Nuclear Data Library
=========================================================

Welcome to the documentation for PapillonNDL ! This is an open source library,
used to read continuous energy ACE files, for neutron cross sections and
secondary distributions. It can find and evaluate cross sections, sample
secondary angle and energy distributions, and provide access to much more data.

Written in C++17 with object orientation, the library provides an easy way to
get the data you need out of ACE files, and fast. A Python API is also generated
using Pybind11, allowing you the same speedy access through Python, without
having to learn a second API! This is a great option for when you want to
make plots of cross sections, or angle/energy distributions.

Currently, only continuous energy neutron data is supported (without photon
production or heating data for now). In the future, we hope to add continuous
energy photon data, as well as thermal neutron scattering data.

The library can be downloaded from GitHub
`here. <https://github.com/HunterBelanger/papillon-ndl>`_ This library is still
young and growing, so if you find a bug or problem please let us know using the
issues feature on GitHub to report the problem. This is by far the most
important way that you can get involved and contribute to the probject ! For
the time being however, pull requests will not be accepted. This will likely
change in the future.


.. toctree::
   :maxdepth: 2

   install/index
   usage/index
   api/index
   license/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
