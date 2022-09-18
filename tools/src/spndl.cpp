#include <docopt.h>

#include <PapillonNDL/absorption.hpp>
#include <PapillonNDL/st_neutron.hpp>
#include <PapillonNDL/elastic_dbrc.hpp>
#include <PapillonNDL/st_coherent_elastic.hpp>
#include <PapillonNDL/st_incoherent_elastic_ace.hpp>
#include <PapillonNDL/st_incoherent_inelastic.hpp>
#include <PapillonNDL/st_thermal_scattering_law.hpp>
#include <cstdint>
#include <iostream>
#include <memory>
#include <ndarray.hpp>
#include <random>

static std::uniform_real_distribution<double> unit(0., 1.);
static std::minstd_rand rng_engine;

double rng() { return unit(rng_engine); }

int reaction(const std::uint32_t mt, const std::uint64_t nsamples,
             const double Ein, const pndl::ACE& ace, NDArray<double>& data) {
  // Get the nuclide
  auto nuclide = pndl::STNeutron(ace);

  // Try to get the reaction in question
  if (nuclide.has_reaction(mt) == false) {
    std::cerr << " Reaction MT " << mt << " is not provided in nuclide.";
    return 1;
  }
  auto& reaction = nuclide.reaction(mt);

  if (Ein <= reaction.threshold()) {
    std::cerr << " Incident energy " << Ein << " MeV is less than threshold of "
              << reaction.threshold() << ".";
    return 1;
  }

  // We have the reaction, now we get the neutron distribution
  auto& distribution = reaction.neutron_distribution();
  if (typeid(distribution).name() == typeid(pndl::Absorption).name()) {
    std::cerr << " Could not find a neutron distribution for MT " << mt << ".";
    return 0;
  }

  // Get all samples
  for (std::uint64_t n = 0; n < nsamples; n++) {
    auto ae = distribution.sample_angle_energy(Ein, rng);

    data(0, n) = ae.energy;
    data(1, n) = ae.cosine_angle;
  }

  return 0;
}

enum class ElasticMode { SVT, DBRC };
int elastic(const ElasticMode mode, const std::uint64_t nsamples,
            const double Ein, const double T, const pndl::ACE& ace,
            NDArray<double>& data) {
  // Get the nuclide
  auto nuclide = pndl::STNeutron(ace);

  nuclide.elastic().set_temperature(T);
  if (mode == ElasticMode::DBRC) {
    nuclide.elastic().set_use_tar(false);
    nuclide.elastic().set_elastic_doppler_broadener(
        std::make_shared<pndl::ElasticDBRC>(nuclide.elastic_xs()));
  }

  // Get all samples
  for (std::uint64_t n = 0; n < nsamples; n++) {
    auto ae = nuclide.elastic().sample_angle_energy(Ein, rng);

    data(0, n) = ae.energy;
    data(1, n) = ae.cosine_angle;
  }

  return 0;
}

int coherent_elastic(const std::uint64_t nsamples, const double Ein,
                     const pndl::ACE& acefile, NDArray<double>& data) {
  pndl::STCoherentElastic distribution(acefile);
  if (distribution.bragg_edges().size() == 0) {
    std::cerr << " TSL does not have coherent elastic scattering.";
    return 1;
  }

  // const auto& distribution = tsl.coherent_elastic();

  // Get all samples
  for (std::uint64_t n = 0; n < nsamples; n++) {
    auto ae = distribution.sample_angle_energy(Ein, rng);

    data(0, n) = ae.energy;
    data(1, n) = ae.cosine_angle;
  }

  return 0;
}

int incoherent_elastic(const std::uint64_t nsamples, const double Ein,
                       const pndl::ACE& acefile, NDArray<double>& data) {
  pndl::STIncoherentElasticACE distribution(acefile);
  if (distribution.cosines().size() == 0) {
    std::cerr << " TSL does not have incoherent elastic scattering.";
    return 1;
  }

  // Get all samples
  for (std::uint64_t n = 0; n < nsamples; n++) {
    auto ae = distribution.sample_angle_energy(Ein, rng);

    data(0, n) = ae.energy;
    data(1, n) = ae.cosine_angle;
  }

  return 0;
}

int incoherent_inelastic(const std::uint64_t nsamples, const double Ein,
                         const pndl::ACE& acefile, NDArray<double>& data) {
  pndl::STIncoherentInelastic distribution(acefile);

  // Get all samples
  for (std::uint64_t n = 0; n < nsamples; n++) {
    auto ae = distribution.sample_angle_energy(Ein, rng);

    data(0, n) = ae.energy;
    data(1, n) = ae.cosine_angle;
  }

  return 0;
}

const std::string HELP_STR =
    "Usage:\n"
    "  spndl reaction <mt> <acefile> <nsamples> <energy> <npyfile>\n"
    "  spndl elastic-svt <acefile> <nsamples> <energy> <T> <npyfile>\n"
    "  spndl elastic-dbrc <acefile> <nsamples> <energy> <T> <npyfile>\n"
    "  spndl coherent-elastic <acefile> <nsamples> <energy> <npyfile>\n"
    "  spndl incoherent-elastic <acefile> <nsamples> <energy> <npyfile>\n"
    "  spndl incoherent-inelastic <acefile> <nsamples> <energy> <npyfile>\n\n"

    "Options:\n"
    "  <mt>        MT identifier of reaction to sample\n"
    "  <acefile>   Name of the ACE file containing data\n"
    "  <nsamples>  Number of samples to perform\n"
    "  <energy>    Incident energy (MeV) at which to take samples\n"
    "  <T>         Temperature of nuclide in Kelvin\n"
    "  <npyfile>   Name of the NPY file in which to write data\n";

enum class RunMode {
  REACTION,
  ELASTIC_SVT,
  ELASTIC_DBRC,
  COHERENT_ELASTIC,
  INCOHERENT_ELASTIC,
  INCOHERENT_INELASTIC
};

int main(int argc, char** argv) {
  // Initialize docopt
  std::map<std::string, docopt::value> args =
      docopt::docopt(HELP_STR, {argv + 1, argv + argc}, false);

  std::uint32_t MT = 0;
  std::uint64_t NSAMPLES = 0;
  double E = 0.;
  double T = 0.;
  std::string ace_file = "";
  std::string npy_file = "";
  RunMode mode = RunMode::REACTION;

  // Get common variables
  NSAMPLES = std::stoul(args["<nsamples>"].asString());
  E = std::stod(args["<energy>"].asString());
  ace_file = args["<acefile>"].asString();
  npy_file = args["<npyfile>"].asString();

  // Get type of reaction to sample
  if (args["reaction"].asBool()) {
    mode = RunMode::REACTION;
    MT = std::stoul(args["<mt>"].asString());
  } else if (args["elastic-svt"].asBool()) {
    mode = RunMode::ELASTIC_SVT;
    T = std::stod(args["<T>"].asString());
  } else if (args["elastic-dbrc"].asBool()) {
    mode = RunMode::ELASTIC_DBRC;
    T = std::stod(args["<T>"].asString());
  } else if (args["coherent-elastic"].asBool()) {
    mode = RunMode::COHERENT_ELASTIC;
  } else if (args["incoherent-elastic"].asBool()) {
    mode = RunMode::INCOHERENT_ELASTIC;
  } else {
    mode = RunMode::INCOHERENT_INELASTIC;
  }

  // Initialize data array
  NDArray<double> data({2, NSAMPLES});
  data.fill(0.);

  // Open ACE File
  pndl::ACE ace(ace_file);

  // Write run info
  std::cout << "==============================================================="
               "======\n";
  std::cout << " PapillonNDL Sampler\n";
  std::cout << " Written by Hunter Belanger\n";
  std::cout << "---------------------------------------------------------------"
               "------\n";
  std::cout << " Input File: " << ace_file << "\n";
  std::cout << " Run Mode: ";
  switch (mode) {
    case RunMode::REACTION:
      std::cout << "Reaction\n";
      break;
    case RunMode::ELASTIC_SVT:
      std::cout << "Elastic SVT\n";
      break;
    case RunMode::ELASTIC_DBRC:
      std::cout << "Elastic DBRC\n";
      break;
    case RunMode::COHERENT_ELASTIC:
      std::cout << "Coherent Elastic\n";
      break;
    case RunMode::INCOHERENT_ELASTIC:
      std::cout << "Incoherent Elastic\n";
      break;
    case RunMode::INCOHERENT_INELASTIC:
      std::cout << "Incoherent Inelastic\n";
      break;
  }
  if (mode == RunMode::REACTION) {
    std::cout << " MT: " << MT << "\n";
  }
  std::cout << " NSAMPLES: " << NSAMPLES << "\n";
  if (mode == RunMode::ELASTIC_SVT || mode == RunMode::ELASTIC_DBRC) {
    std::cout << " Temperature: " << T << " Kelvin\n";
  }
  std::cout << " Energy: " << E << " MeV\n\n";

  if (mode != RunMode::REACTION && E > 4.) {
    std::cerr << "\n WARNING: Sampling a thermal scattering law with an energy "
                 "greater than 4 eV.\n";
    std::cerr << "            Results might not be reliable.\n";
  }

  if (ace.temperature() > 1. && mode == RunMode::ELASTIC_DBRC) {
    std::cerr << "\n WARNING: Sampling Elastic DBRC without 0K elastic xs.\n";
  }

  if (mode == RunMode::ELASTIC_SVT &&
      E >= (400. * (8.617333262E-5) * T * 1.E-6)) {
    std::cerr << "\n WARNING: Asked for SVT, but E > 400 kT. Asymptotic "
                 "approximation will be used.\n";
  }

  // Get Samples
  int result = 0;
  try {
    switch (mode) {
      case RunMode::REACTION:
        result = reaction(MT, NSAMPLES, E, ace, data);
        break;
      case RunMode::ELASTIC_SVT:
        result = elastic(ElasticMode::SVT, NSAMPLES, E, T, ace, data);
        break;
      case RunMode::ELASTIC_DBRC:
        result = elastic(ElasticMode::DBRC, NSAMPLES, E, T, ace, data);
        break;
      case RunMode::COHERENT_ELASTIC:
        result = coherent_elastic(NSAMPLES, E, ace, data);
        break;
      case RunMode::INCOHERENT_ELASTIC:
        result = incoherent_elastic(NSAMPLES, E, ace, data);
        break;
      case RunMode::INCOHERENT_INELASTIC:
        result = incoherent_inelastic(NSAMPLES, E, ace, data);
        break;
    }
  } catch (pndl::PNDLException& error) {
    std::cout << "\n\n !!! ERROR !!!\n\n";
    std::cout
        << " The following problem occured when trying to obtain samples:\n";
    std::cout << error.what() << "\n";
    std::cout << "============================================================="
                 "========\n";
    return 1;
  }

  if (result != 0) {
    std::cout << "\n\n !!! ERROR !!!\n\n";
    std::cout << " Could not generate samples.\n";
  } else {
    std::cout << " Sampling suceeded !\n";
    data.save(npy_file);
    std::cout << " NPY File: " << npy_file << "\n";
  }

  std::cout << "==============================================================="
               "======\n";

  return result;
}
