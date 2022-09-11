#include <docopt.h>
#include <ndarray.hpp>

#include <openmc/endf.h>
#include <openmc/hdf5_interface.h>
#include <openmc/nuclide.h>
#include <openmc/thermal.h>
#include <openmc/constants.h>
#include <openmc/secondary_thermal.h>

#include <iostream>
#include <cstdint>
#include <memory>

using Nuclide = openmc::Nuclide;
using ReactionProduct = openmc::ReactionProduct;
using AngleEnergy = openmc::AngleEnergy;

constexpr std::uint32_t MT_MAX {901};
std::uint64_t seed = 383649624;

std::unique_ptr<openmc::Function1D> elastic_xs = nullptr;
std::unique_ptr<openmc::Function1D> inelastic_xs = nullptr;

std::unique_ptr<Nuclide> get_nuclide(const hid_t h5file) {
  std::unique_ptr<Nuclide> nuclide = nullptr;
  
  // There should only be 1 group in the file. Get the group
  // and then pass it into the Nuclide constructor.
  auto groups = openmc::group_names(h5file);

  hid_t grp_id = openmc::open_group(h5file, groups.at(0));
  nuclide = std::make_unique<Nuclide>(grp_id, std::vector<double>(1, 293.6));
  openmc::close_group(grp_id);

  return nuclide;
}

std::unique_ptr<AngleEnergy> get_coherent_elastic(const hid_t h5file, const std::string& tmpgroup) {
  std::unique_ptr<AngleEnergy> ae = nullptr;

  // There should only be 1 group in the file. Get the group
  // and then pass it into the Nuclide constructor.
  auto groups = openmc::group_names(h5file);
  hid_t grp_id = openmc::open_group(h5file, groups.at(0));
  
  if (openmc::object_exists(grp_id, tmpgroup.data()) == false) {
    std::cerr << " No temperature group '" << tmpgroup << "' exists in file."; 
    openmc::close_group(grp_id);
    return ae;
  }

  hid_t temp_grp_id = openmc::open_group(grp_id, tmpgroup);

  if (openmc::object_exists(temp_grp_id, "elastic")) {
    hid_t elastic_group = openmc::open_group(temp_grp_id, "elastic"); 
    elastic_xs = openmc::read_function(elastic_group, "xs");
    hid_t dgroup = openmc::open_group(elastic_group, "distribution");
    std::string temp;
    openmc::read_attribute(dgroup, "type", temp);
    if (temp == "coherent_elastic") {
      auto xs = dynamic_cast<openmc::CoherentElasticXS*>(elastic_xs.get());
      ae = std::make_unique<openmc::CoherentElasticAE>(*xs);
    } else {
      std::cerr << " No Coherent Elastic provided in file."; 
      openmc::close_group(dgroup);
      openmc::close_group(elastic_group);
      openmc::close_group(temp_grp_id);
      openmc::close_group(grp_id);
      return ae;
    }
    openmc::close_group(dgroup);
    openmc::close_group(elastic_group);
  } else {
    std::cerr << " No Coherent Elastic provided in file.";
    openmc::close_group(temp_grp_id);
    openmc::close_group(grp_id);
    return ae;
  }
  
  openmc::close_group(temp_grp_id);
  openmc::close_group(grp_id);
  return ae;
}

std::unique_ptr<AngleEnergy> get_incoherent_elastic(const hid_t h5file, const std::string& tmpgroup) {
  std::unique_ptr<AngleEnergy> ae = nullptr;

  // There should only be 1 group in the file. Get the group
  // and then pass it into the Nuclide constructor.
  auto groups = openmc::group_names(h5file);
  hid_t grp_id = openmc::open_group(h5file, groups.at(0));
  
  if (openmc::object_exists(grp_id, tmpgroup.data()) == false) {
    std::cerr << " No temperature group '" << tmpgroup << "' exists in file."; 
    openmc::close_group(grp_id);
    return ae;
  }

  hid_t temp_grp_id = openmc::open_group(grp_id, tmpgroup);

  if (openmc::object_exists(temp_grp_id, "elastic")) {
    hid_t elastic_group = openmc::open_group(temp_grp_id, "elastic"); 
    elastic_xs = openmc::read_function(elastic_group, "xs");
    hid_t dgroup = openmc::open_group(elastic_group, "distribution");
    std::string temp;
    openmc::read_attribute(dgroup, "type", temp);
    if (temp == "incoherent_elastic") {
      ae = std::make_unique<openmc::IncoherentElasticAE>(dgroup); 
    } else if (temp == "incoherent_elastic_discrete") {
      auto xs = dynamic_cast<openmc::Tabulated1D*>(elastic_xs.get());
      ae = std::make_unique<openmc::IncoherentElasticAEDiscrete>(dgroup, xs->x());
    } else {
      std::cerr << " No Incoherent Elastic provided in file."; 
      openmc::close_group(dgroup);
      openmc::close_group(elastic_group);
      openmc::close_group(temp_grp_id);
      openmc::close_group(grp_id);
      return ae;
    }
    openmc::close_group(dgroup);
    openmc::close_group(elastic_group);
  } else {
    std::cerr << " No Incoherent Elastic provided in file.";
    openmc::close_group(temp_grp_id);
    openmc::close_group(grp_id);
    return ae;
  }
  
  openmc::close_group(temp_grp_id);
  openmc::close_group(grp_id);
  return ae;
}

std::unique_ptr<AngleEnergy> get_incoherent_inelastic(const hid_t h5file, const std::string& tmpgroup) {
  std::unique_ptr<AngleEnergy> ae = nullptr;

  // There should only be 1 group in the file. Get the group
  // and then pass it into the Nuclide constructor.
  auto groups = openmc::group_names(h5file);
  hid_t grp_id = openmc::open_group(h5file, groups.at(0));
  
  if (openmc::object_exists(grp_id, tmpgroup.data()) == false) {
    std::cerr << " No temperature group '" << tmpgroup << "' exists in file."; 
    openmc::close_group(grp_id);
    return ae;
  }

  hid_t temp_grp_id = openmc::open_group(grp_id, tmpgroup);

  if (openmc::object_exists(temp_grp_id, "inelastic")) {
    hid_t inelastic_group = openmc::open_group(temp_grp_id, "inelastic"); 
    inelastic_xs = openmc::read_function(inelastic_group, "xs");
    hid_t dgroup = openmc::open_group(inelastic_group, "distribution");
    std::string temp;
    openmc::read_attribute(dgroup, "type", temp);
    if (temp == "incoherent_inelastic") {
      ae = std::make_unique<openmc::IncoherentInelasticAE>(dgroup);
    } else if (temp == "incoherent_inelastic_discrete") {
      auto xs = dynamic_cast<openmc::Tabulated1D*>(inelastic_xs.get());
      ae = std::make_unique<openmc::IncoherentInelasticAEDiscrete>(dgroup, xs->x());
    } else {
      std::cerr << " No Incoherent Inelastic provided in file."; 
      openmc::close_group(dgroup);
      openmc::close_group(inelastic_group);
      openmc::close_group(temp_grp_id);
      openmc::close_group(grp_id);
      return ae;
    }
    openmc::close_group(dgroup);
    openmc::close_group(inelastic_group);
  } else {
    std::cerr << " No Incoherent Inelastic provided in file.";
    openmc::close_group(temp_grp_id);
    openmc::close_group(grp_id);
    return ae;
  }
  
  openmc::close_group(temp_grp_id);
  openmc::close_group(grp_id);
  return ae;
}

int reaction(const std::uint32_t mt, const std::uint64_t nsamples,
             const double Ein, const hid_t h5file, NDArray<double>& data) {
  // Get the nuclide
  auto nuclide = get_nuclide(h5file);
  const double AWR = nuclide->awr_;

  // Try to get the reaction in question
  int reaction_index = nuclide->reaction_index_[mt];
  if (reaction_index == openmc::C_NONE) {
    std::cerr << " Reaction MT " << mt << " is not provided in nuclide."; 
    return 1;
  }
  auto& reaction = nuclide->reactions_[reaction_index];
  const bool CM = reaction->scatter_in_cm_;

  if (Ein <= reaction->xs_[0].threshold) {
    std::cerr << " Incident energy " << Ein << " eV is less than threshold of " << reaction->xs_[0].threshold << "."; 
    return 1;
  }

  // We have the reaction, now we get the neutron distribution
  const ReactionProduct* distribution = nullptr;
  for (const auto& product : reaction->products_) {
    if (product.particle_ == openmc::ParticleType::neutron) {
      distribution = &product;
      break;
    }
  }
  if (distribution == nullptr) {
    std::cerr << " Could not find a neutron distribution for MT " << mt << "."; 
    return 0;
  }

  // Get all samples
  for (std::uint64_t n = 0; n < nsamples; n++) {
    double Eout = 0.;
    double mu = -2.; 

    distribution->sample(Ein, Eout, mu, &seed);

    if (CM) {
      double E_cm = Eout; 
      Eout = E_cm + (Ein + 2.0 * mu * (AWR + 1.0) * std::sqrt(Ein * E_cm)) /
                         ((AWR + 1.0) * (AWR + 1.0));
      mu = mu * std::sqrt(E_cm / Eout) + 1.0 / (AWR + 1.0) * std::sqrt(Ein / Eout);
    }

    data(0,n) = Eout / 1.E6; // Convert from eV to MeV for comparison;
    data(1,n) = mu;
  }

  return 0;
}

int coherent_elastic(const std::uint64_t nsamples, const double Ein,
                     const hid_t h5file, const std::string& tmpgroup, NDArray<double>& data) {
  auto distribution = get_coherent_elastic(h5file, tmpgroup);

  if (distribution == nullptr) {
    return 1; 
  }

  // Get all samples
  for (std::uint64_t n = 0; n < nsamples; n++) {
    double Eout = 0.;
    double mu = -2.; 

    distribution->sample(Ein, Eout, mu, &seed);

    data(0,n) = Eout / 1.E6; // Convert from eV to MeV for comparison;
    data(1,n) = mu;
  }

  return 0;
}

int incoherent_elastic(const std::uint64_t nsamples, const double Ein,
                     const hid_t h5file, const std::string& tmpgroup, NDArray<double>& data) {
  auto distribution = get_incoherent_elastic(h5file, tmpgroup);

  if (distribution == nullptr) {
    return 1; 
  }

  // Get all samples
  for (std::uint64_t n = 0; n < nsamples; n++) {
    double Eout = 0.;
    double mu = -2.; 

    distribution->sample(Ein, Eout, mu, &seed);

    data(0,n) = Eout / 1.E6; // Convert from eV to MeV for comparison;
    data(1,n) = mu;
  }

  return 0;
}

int incoherent_inelastic(const std::uint64_t nsamples, const double Ein,
                     const hid_t h5file, const std::string& tmpgroup, NDArray<double>& data) {
  auto distribution = get_incoherent_inelastic(h5file, tmpgroup);

  if (distribution == nullptr) {
    return 1; 
  }

  // Get all samples
  for (std::uint64_t n = 0; n < nsamples; n++) {
    double Eout = 0.;
    double mu = -2.; 

    distribution->sample(Ein, Eout, mu, &seed);

    data(0,n) = Eout / 1.E6; // Convert from eV to MeV for comparison;
    data(1,n) = mu;
  }

  return 0;
}

const std::string HELP_STR = 
  "Usage:\n"
  "  sopenmc reaction <mt> <h5file> <nsamples> <energy> <npyfile>\n"
  "  sopenmc coherent-elastic <h5file> [--temp <tmpgroup>] <nsamples> <energy> <npyfile>\n"
  "  sopenmc incoherent-elastic <h5file> [--temp <tmpgroup>] <nsamples> <energy> <npyfile>\n"
  "  sopenmc incoherent-inelastic <h5file> [--temp <tmpgroup>] <nsamples> <energy> <npyfile>\n\n"
  
  "Options:\n"
  "  <mt>        MT identifier of reaction to sample\n"
  "  <h5file>    Name of the HDF5 file containing data\n"
  "  <tmpgroup>  Temperature group name for HDF5 file (default '294K')\n"
  "  <nsamples>  Number of samples to perform\n"
  "  <energy>    Incident energy (MeV) at which to take samples\n"
  "  <npyfile>   Name of the NPY file in which to write data\n";

enum class RunMode {REACTION, COHERENT_ELASTIC, INCOHERENT_ELASTIC, INCOHERENT_INELASTIC};

int main(int argc, char** argv) {

  // Initialize docopt
  std::map<std::string, docopt::value> args =
    docopt::docopt(HELP_STR, {argv + 1, argv + argc}, false);
  
  std::uint32_t MT = 0;
  std::uint64_t NSAMPLES = 0;
  double E = 0.;
  std::string hdf5_file = "";
  std::string npy_file = "";
  std::string temp_group = "294K";
  RunMode mode = RunMode::REACTION;
  
  // Get common variables  
  NSAMPLES = std::stoul(args["<nsamples>"].asString()); 
  E = std::stod(args["<energy>"].asString()) * 1.E6; // Convert from MeV to eV
  hdf5_file = args["<h5file>"].asString();
  npy_file = args["<npyfile>"].asString();
  if (args["--temp"].asBool()) {
    temp_group = args["<tmpgroup>"].asString();
  }
  
  // Get type of reaction to sample
  if (args["reaction"].asBool()) {
    mode = RunMode::REACTION;
    MT = std::stoul(args["<mt>"].asString());
  } else if (args["coherent-elastic"].asBool()) {
    mode = RunMode::COHERENT_ELASTIC; 
  } else if (args["incoherent-elastic"].asBool()) {
    mode = RunMode::INCOHERENT_ELASTIC; 
  } else {
    mode = RunMode::INCOHERENT_INELASTIC; 
  } 
  
  // Initialize data array
  NDArray<double> data({2,NSAMPLES});
  data.fill(0.);

  // Open HDF5 File
  hid_t h5_file_id = openmc::file_open(hdf5_file, 'r', false);

  // Write run info
  std::cout << "=====================================================================\n";
  std::cout << " OpenMC Sampler\n"; 
  std::cout << " Written by Hunter Belanger\n";
  std::cout << "---------------------------------------------------------------------\n";
  std::cout << " Input File: " << hdf5_file << "\n";
  std::cout << " Run Mode: ";
  switch (mode) {
    case RunMode::REACTION:
      std::cout << "Reaction\n"; 
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
  std::cout << " Energy: " << E*1.E-6 << " MeV\n\n";

  int result = 0;
  if (MT > MT_MAX) {
    std::cerr << " MT must be less than " << MT_MAX << "."; 
    result = 1;
  } else if (MT == 2) {
    std::cerr << " Cannot sample MT 2.";
    result = 1;
  }

  if (mode != RunMode::REACTION && E > 4.) {
    std::cerr << "\n WARNING: Sampling a thermal scattering law with an energy greater than 4 eV.\n"; 
    std::cerr << "            Results might not be reliable.\n";
  }

  // Get Samples
  if (result == 0) {
    switch (mode) {
      case RunMode::REACTION:
        result = reaction(MT, NSAMPLES, E, h5_file_id, data); 
        break;
      case RunMode::COHERENT_ELASTIC:
        result = coherent_elastic(NSAMPLES, E, h5_file_id, temp_group, data);
        break;
      case RunMode::INCOHERENT_ELASTIC:
        result = incoherent_elastic(NSAMPLES, E, h5_file_id, temp_group, data);
        break;
      case RunMode::INCOHERENT_INELASTIC:
        result = incoherent_inelastic(NSAMPLES, E, h5_file_id, temp_group, data);
        break;
    }
  }

  if (result != 0) {
    std::cout << "\n\n !!! ERROR !!!\n\n"; 
    std::cout << " Could not generate samples.\n";
  } else {
    std::cout << " Sampling suceeded !\n";
    data.save(npy_file); 
    std::cout << " NPY File: " << npy_file << "\n";
  }

  std::cout << "=====================================================================\n";
  
  // Close HDF5 File
  openmc::file_close(h5_file_id);

  return result;
}
