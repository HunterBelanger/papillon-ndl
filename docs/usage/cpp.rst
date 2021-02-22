.. _usage_cpp:

===
C++
===

-------------------
Reading an ACE File
-------------------

Most use cases of PapillonNDL only require the inclusion of a single header
file: ``PapillonNDL/nuclide.hpp``. This allows you to first read an ACE file,
and then construct a Nuclide object from it (assuming it is an ACE file with
continuous energy neutron data !).

.. code-block:: c++

  #include <PapillonNDL/nuclide.hpp>
  #include <iostream>

  int main() {

    pndl::ACE U235ace = pndl::ACE("U235.300c");
    pndl::Nuclide U235(U235ace);

    // Writes the ZAID to the terimal. For U235, this should
    // be 92235.
    std::cout << " ZAID: " << U235.ZAID() << std::endl;

    // The temperature of is given in units of kelvin !
    std::cout << " Temp: " << U235.temperature() << std::endl;

    return 0;
  }

When compiling C++ which uses PapillonNDL, be sure to link to the library !
If the header files are already in your include path, and the library is in
your linker search path, then this should be as easy as:

.. code-block:: sh

  g++ pndl_test.cpp -lPapillonNDL -o pndltest

-------------------------
Evaluating Cross Sections
-------------------------

Nuclide objects have quick-access functions to allow fast evaluation of the
total, absorption, and elastic scattering cross sections. If we want to find
these cross sections for U235 at 3MeV, then we can do

.. code-block:: c++

  double tot_xs_at_3mev = U235.total_xs(3.);
  double abs_xs_at_3mev = U235.absorption_xs(3.);
  double ela_xs_at_3mev = U235.elastic_xs(3.);

The standard unit of energy in PapillonNDL is MeV. All energies are given in
MeV, and it expects all arguments which are an energy to be in units of MeV.
While method of evaluating the cross sections works, it is not very efficient.
PapillonNDL has implemented the hashing algorithm implemented in MCNP to speed
up cross section searches. To use this feature, you need to search for the
the index of the desired energy in the EnergyGrid of the nuclide, and then
pass that index to the cross section evaluation call.

.. code-block:: c++

  // Finds index in the energy grid for 3 MeV, using 
  // hashing algorithm for speed.
  size_t i = U235.energy_grid_index(3.);

  tot_xs_at_3mev = U235.total_xs(3., i);
  abs_xs_at_3mev = U235.absorption_xs(3., i);
  ela_xs_at_3mev = U235.elastic_xs(3., i);

This method produces the same results, but is faster overall, and should be
used when performance counts.

If you want the cross section for a specific MT reaction, this can be found
as well. First, it is good practice to see if the nuclide has the desired
reaction

.. code-block:: c++

  // Check to see if U235 has the MT=18 (fission) reaction defined.
  bool has_18 = U235.has_reaction(18);

  double fiss_xs_at_3mev = U235.reaction_xs(18, 3., i);

Passing the index ``i`` is optional, but of course is faster. One can
also skip checking to see if the MT is defined. If it isn't, the
``reaction_xs`` method will return 0.

-------------------------------
Sampling Reaction Distributions
-------------------------------

In Monte Carlo simulations, we often need to sample data related to reactions
such as the number of secondary neutrons produces, and their angle-energy
distributions. To do this, get a reference to the desired reaction; Here, we
will look at the (n,2n) reaction (MT=16):

.. code-block:: c++

  // I know that U235 has MT=16, so we don't need to check that
  // it exists, but this should be done in general !
  const pndl::Reaction& U235_n2n = U236.reaction(16);

  double E_min = U235_n2n.threshold();

  // Here, we get the xs at 6MeV, as 3MeV is bellow the threshold
  // for this reaction !
  double n2n_xs_at_3mev = U235_n2n.xs(6., i);

  double Qval = U235_n2n.Q();

  // For MT=16, the yield is always 2, no matter the energy, but
  // some reactions has energy dependent yields.
  double n_out = U235_n2n.yield(6.); 

In the above example, we have been able to get lots of data about the
reaction, such as the Q-value, the minimum energy at which is occurs,
and the reaction channels yield. Before we can sample from the secondary
distributions however, we need a random number generator function, which
produces random doubles on the interval [0,1). We will set one up really
fast to demonstrate how sampling works.

.. code-block:: c++

  #include <random>

  std::minstd_rand rng_eng;
  std::uniform_real_distribution<> U(0.,1.);

  double rng() {
    return U(rng_eng);
  }

This isn't exactly beautiful, but it gets the job done. A random number
generator function must be provided as some of the algorithms to sample
the energy distributions require many random numbers, and it is
impossible to know how many it will need in advance. We can now sample
and outgoing angle and energy in the laboratory frame with

.. code-block:: c++

  pndl::AngleEnergyPacket out = U235_n2n.sample_angle_energy(6., rng);

The cose of the scattering angle is then stored in ``out.angle``, and
the energy is in ``out.energy``.

------------
Fission Data
------------

Often we want to look up lots of particular fission data for isotopes
such as U235. While the fission cross section is contained in the
MT=18 reaction, the best place to get other bits of fission data such
as the number of neutrons per fission, the fission neutron spectrum,
and delayed group info/spectra is the FissionData object in the
nuclide.

.. code-block:: c++

  const pndl::FissionData& fiss_data = U235.fission_data();

  // Total number of fission neutrons for fissions induced by 3 MeV
  // neutrons.
  double nu = fiss_data.nu_total(3.);

  double nu_prmpt = fiss_data.nu_prompt(3.);
  double nu_delyd = fiss_data.nu_delayed(3.);

  // The angle-energy distributions for prompt neutrons can be sampled
  // with the method
  pndl::AngleEnergyPacket fiss_out = fiss_data.sample_prompt_angle_energy(3., rng);

Information for a delayed neutron group is also available in a DelayedGroup class:

.. code-block:: c++

  size_t delayed_grps = fiss_data.ngroups();

  // Delayed groups are indexed starting from 0
  const pndl::DelayedGroup& dg1 = fiss_data.delayed_group(1);

  // The decay constant for the group is given in units of
  // inverse seconds.
  double decay_const = dg1.decay_constant();

  // The probability of a fission neutron being in the given group is a
  // function of the incident energy
  double prob_dg1 = dg1.probability(3.);

It is always assumed that the neutrons born from a delayed group have an
isotropic angular distribution. As such, we only sample the energy from
the delayed group

.. code-block:: c++

  double E_out = dg1.sample_energy(3., rng);

This should be enough of an introduction for most users to start using the
library to get work done, and access continuous energy neutron data. In an
effort to maintain SOLID programming principles, energy and angle distributions
are only ever accesed through virtual interface classes. This is not the case
for the Python bindings however, as these are generated with Pybind11, which
always downcasts pointers the pointers to the true object type. This makes it
possible to see more of the inner workings of the library, and gain access to
specific parts of distributions. If this is what you're into, take a look at
using the Python API. It's just as fast (as it is in C++), but is a gereat
way to plot data, especially cross sections and distributions.