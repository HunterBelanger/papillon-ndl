/*
 * Copyright 2021, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *
 * */
#include <PapillonNDL/ce_neutron.hpp>
#include <PapillonNDL/constant.hpp>
#include <PapillonNDL/multi_region_1d.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/polynomial_1d.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <PapillonNDL/sum_1d.hpp>
#include <PapillonNDL/difference_1d.hpp>
#include <PapillonNDL/uncorrelated.hpp>
#include <memory>
#include <system_error>
#include <vector>

namespace pndl {

CENeutron::CENeutron(const ACE& ace)
    : zaid_(ace.zaid()),
      awr_(ace.awr()),
      temperature_(ace.temperature()),
      fissile_(ace.fissile()),
      energy_grid_(nullptr),
      total_xs_(nullptr),
      disappearance_xs_(nullptr),
      elastic_xs_(nullptr),
      fission_xs_(nullptr),
      photon_production_xs_(nullptr),
      elastic_angle_(nullptr),
      nu_total_(nullptr),
      nu_prompt_(nullptr),
      nu_delayed_(nullptr),
      delayed_groups_(),
      mt_list_(),
      reaction_indices_(),
      reactions_() {
  // Construct energy grid
  energy_grid_ = std::make_shared<EnergyGrid>(ace);

  // Number of energy points
  uint32_t NE = ace.nxs(2);
  total_xs_ =
      std::make_shared<CrossSection>(ace, ace.ESZ() + NE, energy_grid_, false);
  disappearance_xs_ = std::make_shared<CrossSection>(ace, ace.ESZ() + 2 * NE,
                                                     energy_grid_, false);
  elastic_xs_ = std::make_shared<CrossSection>(ace, ace.ESZ() + 3 * NE,
                                               energy_grid_, false);

  // Get photon production XS if present
  if (ace.jxs(11) != 0) {
    photon_production_xs_ =
        std::make_shared<CrossSection>(ace, ace.GPD(), energy_grid_, false);
  } else {
    photon_production_xs_ = std::make_shared<CrossSection>(0., energy_grid_); 
  }

  // Make elastic AngleDistribution
  elastic_angle_ =
      std::make_shared<AngleDistribution>(ace, ace.xss<int>(ace.LAND()));

  // Read all reactions
  uint32_t NMT = ace.nxs(3);
  reaction_indices_.fill(-1);
  int32_t current_reaction_index = 0;
  mt_list_.resize(NMT, 0);
  reactions_.reserve(NMT);
  for (uint32_t indx = 0; indx < NMT; indx++) {
    uint32_t MT = ace.xss<uint32_t>(ace.MTR() + indx);
    mt_list_[indx] = MT;
    reactions_.emplace_back(ace, indx, energy_grid_);
    reaction_indices_[MT] = current_reaction_index;
    current_reaction_index++;
  }

  if (fissile_) {
    read_fission_data(ace);
  } else {
    nu_total_ = std::make_shared<Constant>(0.);
    nu_prompt_ = std::make_shared<Constant>(0.);
    nu_delayed_ = std::make_shared<Constant>(0.);
  }
  
  fission_xs_ = compute_fission_xs(); 
}

CENeutron::CENeutron(const ACE& ace, const CENeutron& nuclide)
    : zaid_(ace.zaid()),
      awr_(ace.awr()),
      temperature_(ace.temperature()),
      fissile_(ace.fissile()),
      energy_grid_(nullptr),
      total_xs_(nullptr),
      disappearance_xs_(nullptr),
      elastic_xs_(nullptr),
      fission_xs_(nullptr),
      photon_production_xs_(nullptr),
      elastic_angle_(nullptr),
      nu_total_(nullptr),
      nu_prompt_(nullptr),
      nu_delayed_(nullptr),
      delayed_groups_(),
      mt_list_(),
      reaction_indices_(),
      reactions_() {
  // Construct energy grid
  energy_grid_ = std::make_shared<EnergyGrid>(ace);

  // Make sure these are the same nuclide !
  if (zaid_ != nuclide.zaid()) {
    std::string mssg =
        "Nuclide::Nuclide: ZAID of ACE doesn't match ZAID of nuclide.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (awr_ != nuclide.awr()) {
    std::string mssg =
        "Nuclide::Nuclide: AWR of ACE doesn't match AWR of nuclide.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Number of energy points
  uint32_t NE = ace.nxs(2);

  total_xs_ =
      std::make_shared<CrossSection>(ace, ace.ESZ() + NE, energy_grid_, false);
  disappearance_xs_ = std::make_shared<CrossSection>(ace, ace.ESZ() + 2 * NE,
                                                     energy_grid_, false);
  elastic_xs_ = std::make_shared<CrossSection>(ace, ace.ESZ() + 3 * NE,
                                               energy_grid_, false);

  // Get photon production XS if present
  if (ace.jxs(11) != 0) {
    photon_production_xs_ =
        std::make_shared<CrossSection>(ace, ace.GPD(), energy_grid_, false);
  }

  // Copy elastic AngleDistribution
  elastic_angle_ = nuclide.elastic_angle_;

  // Read all reactions
  uint32_t NMT = ace.nxs(3);
  mt_list_.resize(NMT, 0);
  reaction_indices_.fill(-1);
  int32_t current_reaction_index = 0;
  reactions_.reserve(NMT);
  for (uint32_t indx = 0; indx < NMT; indx++) {
    uint32_t MT = ace.xss<uint32_t>(ace.MTR() + indx);
    mt_list_[indx] = MT;

    if (!nuclide.has_reaction(MT)) {
      std::string mssg = "Nuclide::Nuclide: MT = " + std::to_string(MT) +
                         " is present in ACE, but not in nuclide.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }

    reactions_.emplace_back(ace, indx, energy_grid_, nuclide.reaction(MT));
    reaction_indices_[MT] = current_reaction_index;
    current_reaction_index++;
  }

  // Copy fission data from other nuclide
  nu_total_ = nuclide.nu_total_;
  nu_prompt_ = nuclide.nu_prompt_;
  nu_delayed_ = nuclide.nu_delayed_;
  delayed_groups_ = nuclide.delayed_groups_;

  fission_xs_ = compute_fission_xs();
}

void CENeutron::read_fission_data(const ACE& ace) {
  // These will temporarily hold the nu functions untill all are made
  // and we can make the Sum1D or Difference1D instances from them.
  std::shared_ptr<Function1D> total(nullptr);
  std::shared_ptr<Function1D> prompt(nullptr);
  std::shared_ptr<Function1D> delayed(nullptr);

  // If prompt and or total neutrons are given
  if (ace.jxs(1) > 0) {
    if (ace.xss(ace.NU()) > 0.) {
      // Either prompt or total given, but not both
      if (ace.jxs(23) > 0) {
        // Prompt is provided, as delayed is present
        prompt = read_nu(ace, ace.NU());
      } else {
        // No delayed, so this must be total
        total = read_nu(ace, ace.NU());
      }
    } else {
      // Both prompt and total given
      uint32_t KNU_prmpt = ace.NU() + 1;
      uint32_t KNU_tot = ace.NU() + std::abs(ace.xss<int32_t>(ace.NU())) + 1;

      total = read_nu(ace, KNU_tot);
      prompt = read_nu(ace, KNU_prmpt);
    }
  }

  // Read delayed nu if given
  if (ace.DNU() > 0) {
    delayed = read_nu(ace, ace.DNU());
  }

  // First we make nu_total_
  if (total) {
    nu_total_ = total; 
  } else if (prompt && delayed) {
    nu_total_ = std::make_shared<Sum1D>(prompt, delayed); 
  } else if (delayed) {
    nu_total_ = delayed; 
  } else {
    // Should never get here, but just in case 
    nu_total_ = std::make_shared<Constant>(0.);
  }

  // Now we make nu_prompt_
  if (prompt) {
    nu_prompt_ = prompt; 
  } else {
    nu_prompt_ = nu_total_; 
  }

  // And finally, delayed
  if (delayed) {
    nu_delayed_ = delayed; 
  } else if(total && prompt) {
    nu_delayed_ = std::make_shared<Difference1D>(total, prompt); 
  } else {
    nu_delayed_ = std::make_shared<Constant>(0.); 
  }

  // Read all delayed group data
  if (ace.BDD() > 0) {
    uint32_t NGRPS = ace.nxs(7);
    std::size_t g = 1;
    std::size_t i = ace.BDD();
    while (g <= NGRPS) {
      delayed_groups_.push_back(DelayedGroup(ace, i, g));
      uint32_t NR = ace.xss<uint32_t>(i + 1);
      uint32_t NE = ace.xss<uint32_t>(i + 2 + 2 * NR);
      i += 3 + 2 * (NR + NE);
      g++;
    }
  }
}

std::shared_ptr<Function1D> CENeutron::read_nu(const ACE& ace, std::size_t i) {
  uint32_t LNU = ace.xss<uint32_t>(i);

  if (LNU == 1) {  // Polynomial
    return read_polynomial_nu(ace, ++i);
  } else {  // Tabular
    return read_tabular_nu(ace, ++i);
  }
}

std::shared_ptr<Function1D> CENeutron::read_polynomial_nu(const ACE& ace,
                                                          std::size_t i) {
  uint32_t NC = ace.xss<uint32_t>(i);
  std::vector<double> coeffs = ace.xss(i + 1, NC);
  return std::make_shared<Polynomial1D>(coeffs);
}

std::shared_ptr<Function1D> CENeutron::read_tabular_nu(const ACE& ace,
                                                       std::size_t i) {
  uint32_t NR = ace.xss<uint32_t>(i);
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<double> energy = ace.xss(i + 2 + 2 * NR, NE);
  std::vector<double> y = ace.xss(i + 2 + 2 * NR + NE, NE);

  if (NR == 0 || NR == 1) {
    Interpolation interp = Interpolation::LinLin;
    if (NR == 1) interp = ace.xss<Interpolation>(i + 2);

    return std::make_shared<Region1D>(energy, y, interp);
  } else {
    std::vector<uint32_t> breaks = ace.xss<uint32_t>(i + 1, NR);
    std::vector<Interpolation> interps = ace.xss<Interpolation>(i + 1 + NR, NR);
    return std::make_shared<MultiRegion1D>(breaks, interps, energy, y);
  }
}

std::shared_ptr<CrossSection> CENeutron::compute_fission_xs() {
  if(!fissile_) {
    return std::make_shared<CrossSection>(0., energy_grid_); 
  }

  if (this->has_reaction(18)) {
    return std::make_shared<CrossSection>(this->reaction(18).xs()); 
  }

  // Life is difficult. We need to sum the products.
  std::size_t lowest_index = energy_grid_->size(); 
  if (this->has_reaction(19)) {
    if (this->reaction(19).xs().index() < lowest_index) {
      lowest_index = this->reaction(19).xs().index(); 
    } 
  }

  if (this->has_reaction(20)) {
    if (this->reaction(20).xs().index() < lowest_index) {
      lowest_index = this->reaction(20).xs().index(); 
    } 
  }

  if (this->has_reaction(21)) {
    if (this->reaction(21).xs().index() < lowest_index) {
      lowest_index = this->reaction(21).xs().index(); 
    } 
  }

  if (this->has_reaction(38)) {
    if (this->reaction(38).xs().index() < lowest_index) {
      lowest_index = this->reaction(38).xs().index(); 
    } 
  }

  // We now know what energy point to start at. We can initialize the vector.
  std::vector<double> fiss_xs(energy_grid_->size() - lowest_index, 0.);

  for (std::size_t i = lowest_index; i < energy_grid_->size(); i++) {
    if (this->has_reaction(19)) {
      fiss_xs[i - lowest_index] += this->reaction(19).xs()[i]; 
    } 

    if (this->has_reaction(20)) {
      fiss_xs[i - lowest_index] += this->reaction(20).xs()[i]; 
    }

    if (this->has_reaction(21)) {
      fiss_xs[i - lowest_index] += this->reaction(21).xs()[i]; 
    }

    if (this->has_reaction(38)) {
      fiss_xs[i - lowest_index] += this->reaction(38).xs()[i]; 
    }
  }

  return std::make_shared<CrossSection>(fiss_xs, energy_grid_, lowest_index);
}

}  // namespace pndl
