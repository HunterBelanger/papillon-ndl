/*
 * Papillon Nuclear Data Library
 * Copyright 2021, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Papillon Nuclear Data Library (PapillonNDL).
 *
 * PapillonNDL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PapillonNDL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PapillonNDL. If not, see <https://www.gnu.org/licenses/>.
 *
 * */
#include <PapillonNDL/ce_neutron.hpp>
#include <PapillonNDL/elastic_svt.hpp>
#include <PapillonNDL/pndl_exception.hpp>

namespace pndl {

CENeutron<CrossSection>::CENeutron(const ACE& ace)
    : CENeutronBase(ace),
      temperature_(ace.temperature()),
      energy_grid_(nullptr),
      total_xs_(nullptr),
      disappearance_xs_(nullptr),
      elastic_xs_(nullptr),
      heating_number_(nullptr),
      fission_xs_(nullptr),
      photon_production_xs_(nullptr),
      elastic_distribution_(nullptr),
      reactions_(),
      urr_ptables_(nullptr) {
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
  heating_number_ = std::make_shared<CrossSection>(ace, ace.ESZ() + 4 * NE,
                                                   energy_grid_, false);

  // Get photon production XS if present
  if (ace.jxs(11) != 0) {
    photon_production_xs_ =
        std::make_shared<CrossSection>(ace, ace.GPD(), energy_grid_, false);
  } else {
    photon_production_xs_ = std::make_shared<CrossSection>(0., energy_grid_);
  }

  // Read all reactions
  uint32_t NMT = ace.nxs(3);
  reaction_indices_.fill(-1);
  int32_t current_reaction_index = 0;
  mt_list_.resize(NMT, 0);
  reactions_.reserve(NMT);
  for (uint32_t indx = 0; indx < mt_list_.size(); indx++) {
    uint32_t MT = ace.xss<uint32_t>(ace.MTR() + indx);
    mt_list_[indx] = MT;
    reactions_.emplace_back(ace, indx, energy_grid_);
    reaction_indices_[MT] = current_reaction_index;
    current_reaction_index++;
  }

  fission_xs_ = compute_fission_xs();

  // Make the PTables. Grab reference to MT 102
  if (this->has_reaction(102) == false) {
    std::string mssg =
        "Nuclide does not have a radiative capture cross section.";
    throw PNDLException(mssg);
  }
  std::shared_ptr<CrossSection> capture_xs_ =
      std::make_shared<CrossSection>(this->reaction(102).xs());
  try {
    urr_ptables_ = std::make_shared<URRPTables>(ace, *elastic_xs_, *capture_xs_,
                                                *fission_xs_, *heating_number_,
                                                reactions_);
  } catch (PNDLException& error) {
    std::string mssg = "Could not construct URRPTables for nuclide data.";
    error.add_to_exception(mssg);
    throw error;
  }

  // Make elastic scatter distribution
  try {
    elastic_distribution_ = std::make_shared<Elastic>(
        std::make_shared<ElasticSVT>(), elastic_angle_, awr_, temperature_);
  } catch (PNDLException& err) {
    std::string mssg = "Could not create Elastic AngleEnergy distribution.";
    err.add_to_exception(mssg);
    throw err;
  }
}

CENeutron<CrossSection>::CENeutron(const ACE& ace, const CENeutron& nuclide)
    : CENeutronBase(nuclide),
      temperature_(ace.temperature()),
      energy_grid_(nullptr),
      total_xs_(nullptr),
      disappearance_xs_(nullptr),
      elastic_xs_(nullptr),
      heating_number_(nullptr),
      fission_xs_(nullptr),
      photon_production_xs_(nullptr),
      elastic_distribution_(nullptr),
      reactions_(),
      urr_ptables_(nullptr) {
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
  heating_number_ = std::make_shared<CrossSection>(ace, ace.ESZ() + 4 * NE,
                                                   energy_grid_, false);

  // Get photon production XS if present
  if (ace.jxs(11) != 0) {
    photon_production_xs_ =
        std::make_shared<CrossSection>(ace, ace.GPD(), energy_grid_, false);
  } else {
    photon_production_xs_ = std::make_shared<CrossSection>(0., energy_grid_);
  }

  // Read all reactions
  uint32_t NMT = ace.nxs(3);
  reaction_indices_.fill(-1);
  int32_t current_reaction_index = 0;
  mt_list_.resize(NMT, 0);
  reactions_.reserve(NMT);
  for (uint32_t indx = 0; indx < mt_list_.size(); indx++) {
    uint32_t MT = ace.xss<uint32_t>(ace.MTR() + indx);
    mt_list_[indx] = MT;
    reactions_.emplace_back(ace, indx, energy_grid_,
                            nuclide.reactions_[nuclide.reaction_indices_[MT]]);
    reaction_indices_[MT] = current_reaction_index;
    current_reaction_index++;
  }

  fission_xs_ = compute_fission_xs();

  // Make the PTables. Grab reference to MT 102
  if (this->has_reaction(102) == false) {
    std::string mssg =
        "Nuclide does not have a radiative capture cross section.";
    throw PNDLException(mssg);
  }
  std::shared_ptr<CrossSection> capture_xs_ =
      std::make_shared<CrossSection>(this->reaction(102).xs());
  try {
    urr_ptables_ = std::make_shared<URRPTables>(ace, *elastic_xs_, *capture_xs_,
                                                *fission_xs_, *heating_number_,
                                                reactions_);
  } catch (PNDLException& error) {
    std::string mssg = "Could not construct URRPTables for nuclide data.";
    error.add_to_exception(mssg);
    throw error;
  }

  // Make elastic scatter distribution
  try {
    elastic_distribution_ = nuclide.elastic_distribution_->clone();
    elastic_distribution_->set_temperature(temperature_);
  } catch (PNDLException& err) {
    std::string mssg = "Could not create Elastic AngleEnergy distribution.";
    err.add_to_exception(mssg);
    throw err;
  }
}

std::shared_ptr<CrossSection> CENeutron<CrossSection>::compute_fission_xs() {
  if (!fissile_) {
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
