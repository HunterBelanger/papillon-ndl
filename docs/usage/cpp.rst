.. _usage_cpp:

===
C++
===

-------------------
Reading an ACE File
-------------------

Most use cases of PapillonNDL only require the inclusion of a single header
file: ``PapillonNDL/st_neutron.hpp``. This allows you to first read an ACE file,
and then construct a STNeutron object from it (assuming it is an ACE file with
continuous energy neutron data!).

.. code-block:: c++

  #include <PapillonNDL/st_neutron.hpp>
  #include <iostream>

  int main() {
    // Read ACE file with nuclear data
    pndl::ACE U235ace = pndl::ACE("U235.300c");

    // Construct an STNeutron instance from the ACE file.
    // And STNeutron contains all continuous energy neutron data
    // for a single nuclide, and at a single temperature.
    pndl::STNeutron U235(U235ace);

    // Writes the ZAID to the terminal. For U235, this should
    // be 92235.
    std::cout << " ZAID: " << U235.zaid() << std::endl;

    // The temperature is given in units of kelvin !
    std::cout << " Temp: " << U235.temperature() << std::endl;

    return 0;
  }

When compiling C++ which uses PapillonNDL, be sure to link to the library!
If the header files are already in your include path, and the library is in
your linker search path, then this should be as easy as:

.. code-block:: sh

  g++ --std=c++20 pndl_test.cpp -lPapillonNDL -o pndltest

------------------------------
Loading a Nuclear Data Library
------------------------------

If you already have an ACE library on your system, from either MCNP or Serpent,
you can start using this data easily with the MCNPLibrary and SerpentLibrary
helper classes. Including the relative header file is all is necessary. As an
example, if you want to use an MCNP xsdir file, then first make sure that the
``PapillonNDL/mcnp_library.hpp`` header file included in the executable.

.. code-block:: c++

  pndl::MCNPLibrary lib("/home/hunter/Documents/nuclear_data/lib80x/ace/lib80x.xsdir");

From a library, you can directly load a nuclide using the element symbol, and
the atomic mass number of the isotope. You must also specify the desired
temperature for the data, and a temperature tolerance (which is 1 Kelvin by
default). If you don’t know what temperatures are provided in the library, you
can get a list of the provided data for the symbol:

.. code-block:: c++

   const std::vector<double>& U238_temps = lib.temperatures("U238");

Let's grab the evaluation for U238 at 293.6 Kelvin:

.. code-block:: c++

   std::shared_ptr<pndl::STNeutron> U238 = lib.load_STNeutron("U238", 293.6);

A library will return the STNeutron in a shared pointer, allowing it to keep
its own copy, so that the data wont be reloaded and dupilcated if the same
nuclide and temperature combination is requested later on.

-------------------------
Evaluating Cross Sections
-------------------------

STNeutron objects have quick-access functions to allow fast evaluation of the
total, absorption, and elastic scattering cross sections. If we want to find
these cross sections for U235 at 3MeV, then we can do

.. code-block:: c++

  double tot_xs_at_3mev = U235.total_xs()(3.);
  double abs_xs_at_3mev = U235.disappearance_xs()(3.);
  double ela_xs_at_3mev = U235.elastic_xs()(3.);

The standard unit of energy in PapillonNDL is MeV. All energies are given in
MeV, and it expects all arguments which are an energy to be in units of MeV.
While this method of evaluating the cross sections works, it is not very efficient.
PapillonNDL has implemented the hashing algorithm implemented in MCNP to speed
up cross section searches. To use this feature, you need to search for the
the index of the desired energy in the EnergyGrid of the nuclide, and then
pass that index to the cross section evaluation call.

.. code-block:: c++

  // Finds index in the energy grid for 3 MeV, using 
  // hashing algorithm for speed.
  size_t i = U235.energy_grid().get_lower_index(3.);

  tot_xs_at_3mev = U235.total_xs()(3., i);
  abs_xs_at_3mev = U235.disappearance_xs()(3., i);
  ela_xs_at_3mev = U235.elastic_xs()(3., i);

This method produces the same results, but is faster overall, and should be
used when performance counts.

If you want the cross section for a specific MT reaction, this can be found
as well. First, it is good practice to see if the nuclide has the desired
reaction

.. code-block:: c++

  // Check to see if U235 has the MT=18 (fission) reaction defined.
  bool has_18 = U235.has_reaction(18);

  double fiss_xs_at_3mev = U235.reaction(18).xs()(3., i);

Passing the index ``i`` is optional, but is of course faster.

-------------------------------
Sampling Reaction Distributions
-------------------------------

In Monte Carlo simulations, we often need to sample data related to reactions
such as the number of secondary neutrons produced, and their angle-energy
distributions. To do this, start by getting a reference to the desired
reaction; here, we will look at the (n,2n) reaction (MT=16):

.. code-block:: c++

  // I know that U235 has MT=16, so we don't need to check that
  // it exists, but this should be done in general !
  const pndl::STReaction& U235_n2n = U235.reaction(16);

  double E_min = U235_n2n.threshold();

  // Here, we get the xs at 6MeV, as 3MeV is below the threshold
  // for this reaction !
  double n2n_xs_at_3mev = U235_n2n.xs()(6., i);

  double Qval = U235_n2n.q();

  // For MT=16, the yield is always 2, no matter the energy, but
  // some reactions have energy dependent yields.
  double n_out = U235_n2n.yield()(6.); 

In the above example, we have been able to get lots of data about the
reaction, such as the Q-value, the minimum energy tabulated for the reaction,
and the reaction channel's yield. Before we can sample from the secondary
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
an outgoing angle and energy in the laboratory frame with

.. code-block:: c++

  pndl::AngleEnergyPacket out = U235_n2n.sample_neutron_angle_energy(6., rng);

The cosine of the scattering angle is then stored in ``out.cosine_angle``,
and the energy is in ``out.energy``.

Absorption reactions which do not emit neutrons have a special type of
distribution which will throw a PNDLException if you try to sample them.

------------------
Elastic Scattering
------------------

PapillonNDL treats elastic scattering differently than other reactions. This is
because many different algorithms can be used for elastic scattering, and the
choice of algorithm can have a large impact on simulation results. The elastic
scattering data is stored in the Elastic class, which also inherits from
AngleEnergy.

By default, PapillonNDL will use the Sample Velocity of Target (SVT) method to
sample elastic scattering. This approximation is also refered to as the
Constant Cross Section (CXS) approximation. While used ubiquitously in Monte
Carlo codes, it is known to give inaccurate results when used for large
nucleids which have resonances at low energies. If desired, you can manally
change the approximation, by giving the Elastic instance a new
ElasticDopplerBroadener. The two possible broadeners are ElasticSVT (the
default), or ElasticDBRC, which applied the Doppler Broadening Rejection
Correction (DBRC). This method requires the 0 Kelvin elastic scattering cross
section, so we will load that, and then apply DBRC to U235.

.. code-block:: c++

   std::shared_ptr<pndl::STNeutron> U235_0K = lib.load_STNeutron("U235", 0.);
   auto dbrc = std::make_shared<pndl::ElasticDBRC>(U235_0K->elastic_xs());
   U235.elastic().set_elastic_doppler_broadener(dbrc);

Another approximation we can change is the use of the Target at Rest (TAR)
approximation. By default, TAR is used for all nuclides when the incident
energy Ein is larger than 400kT, where k is the Boltzmann constant, and T is
the nuclide temperature. It is generally a good idea to use this approximation
when Ein > 400kT, as it does not significantly change results, and speeds up
calculations. Most codes however do not use TAR for Hydrogen-1, as it is
lighter than a neutron, and using this approximation can lead to inaccurate
energy transfers. You can turn off TAR with

.. code-block:: c++

   U235.elastic().set_use_tar(false);

You may prefer to change the threshold at which TAR is applied. If you would
rather only use TAR when Ein > 700kT, then you can use something like the
following:

.. code-block:: c++

   U235.elastic().set_use_tar(true);
   U235.elastic().set_tar_threshold(700.);

------------
Fission Data
------------

Often we want to look up lots of particular fission data for isotopes
such as U235. While the fission cross section is directly stored in the
STNeutron, other pieces of fission data are stored in the Fission class,
contained in the STNeutron instance. Methods of the Fission class allow us to
access other bits of fission data such as the number of neutrons per
fission, the prompt neutron spectrum, and delayed family info/spectra.

.. code-block:: c++

  // Total number of fission neutrons for fissions induced by 3 MeV
  // neutrons.
  double nu = U235.fission().nu_total()(3.);
  double nu_prmpt = U235.fission().nu_prompt()(3.);
  double nu_delyd = U235.fission().nu_delayed()(3.);

The prompt spectrum is also provided here, and can be sampled like a regular
reaction distribution.

.. code-block:: c++
  
  // Sample an angle-energy pair for prompt fission induced at 1.2 eV.
  pndl::AngleEnergyPacket prmpt = U235.fission().prompt_spectrum().sample_angle_energy(1.2E-6, rng);

Information for a delayed neutron family is also available in a DelayedFamily class:

.. code-block:: c++

  size_t delayed_fms = U235.fission().n_delayed_families();

  // Delayed families are indexed starting from 0
  const pndl::DelayedFamily& df1 = U235.fission().delayed_family(1);

  // The decay constant for the family is given in units of
  // inverse seconds.
  double decay_const = df1.decay_constant();

  // The probability of a fission neutron being in the given family is a
  // function of the incident energy
  double prob_df1 = df1.probability()(3.);

It is always assumed that the neutrons born from a delayed family have an
isotropic angular distribution. As such, we only sample the energy from
the delayed family.

.. code-block:: c++

  double E_out = df1.sample_energy(3., rng);

---------------------------
Unresolved Resonance Region
---------------------------

At high energies, it becomes impossible to determine resonanace parameters.
This region is called the Unresolved Resonance Region (URR). To correctly treat
this portion of the energy spectrum, probability tables can be used. This
information is stored in the URRPTables class. If the result of
URRPTables.is_valid() is true, then URR data is provided, and can be used in
transport.

.. code-block:: c++

   if (U235.urr_ptables().is_valid()) {
     std::cout << "URR PTables are valid for U235 !\n";
   }

We can get the minimum and maximum energy for the URR with the following
methods:

.. code-block:: c++

   double URR_Emin = U235.urr_ptables().min_energy();
   double URR_Emax = U235.urr_ptables().max_energy();

If we have a neutron energy which is within the URR region, we sample a random
number xi. This same random value must be used for the given nuclide, no matter
the temperature, until the neutron has undergone a collision. It is used to
sample the cross sections from the probability tables, for the given nuclide
over the given flight.

.. code-block:: c++

   double Ein = 0.01; // 0.01 MeV
   double xi = rng();
   pndl::XSPacket xs;
   
   if (URR_Emin < Ein && Ein < URR_Emax) {
     xs = U235.urr_ptables().evaluate_xs(Ein, xi).value();

     std::cout << "Total XS = " << xs.total << "\n";
     std::cout << "Elasic XS = " << xs.elastic << "\n";
     std::cout << "Inelastic XS = " << xs.inelastic << "\n";
     std::cout << "Absorption XS = " << xs.absorption << "\n";
     std::cout << "Fission XS = " << xs.fission << "\n";
   } 
   

This should be enough of an introduction for most users to start using the
library to get work done, and access continuous energy neutron data. In an
effort to maintain SOLID programming principles, energy and angle distributions
are only ever accesed through virtual interface classes. This is not the case
for the Python bindings however, as these are generated with Pybind11, which
always downcasts objects to the true type. This makes it possible to see more
of the inner workings of the library, and gain access to specific parts of
distributions. If this is what you're into, take a look at using the Python
API. It's just as fast (as it is written in C++), but is a gereat way to plot
data, especially cross sections and distributions.
