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

namespace pndl {

CENeutron<CrossSection>::CENeutron(const ACE& ace)
    : CENeutronBase(ace),
      temperature_(ace.temperature()),
      energy_grid_(ace),
      total_xs_(nullptr),
      disappearance_xs_(nullptr),
      elastic_xs_(nullptr),
      fission_xs_(nullptr),
      photon_production_xs_(nullptr),
      reactions_() {
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
}

CENeutron<CrossSection>::CENeutron(const ACE& ace, const CENeutron& nuclide)
    : CENeutronBase(nuclide),
      temperature_(ace.temperature()),
      energy_grid_(ace),
      total_xs_(nullptr),
      disappearance_xs_(nullptr),
      elastic_xs_(nullptr),
      fission_xs_(nullptr),
      photon_production_xs_(nullptr),
      reactions_() {
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
}

std::shared_ptr<CrossSection> CENeutron<CrossSection>::compute_fission_xs() {
  if (!fissile_) {
    return std::make_shared<CrossSection>(0., energy_grid_);
  }

  if (this->has_reaction(18)) {
    return std::make_shared<CrossSection>(this->reaction(18).xs());
  }

  // Life is difficult. We need to sum the products.
  std::size_t lowest_index = energy_grid_.size();
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
  std::vector<double> fiss_xs(energy_grid_.size() - lowest_index, 0.);

  for (std::size_t i = lowest_index; i < energy_grid_.size(); i++) {
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
