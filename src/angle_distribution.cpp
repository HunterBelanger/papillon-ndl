/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2022, Hunter Belanger
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
#include <PapillonNDL/angle_distribution.hpp>
#include <PapillonNDL/angle_table.hpp>
#include <PapillonNDL/equiprobable_angle_bins.hpp>
#include <PapillonNDL/isotropic.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <cmath>
#include <memory>

namespace pndl {

AngleDistribution::AngleDistribution() : energy_grid_(), laws_() {
  energy_grid_.reserve(2);
  laws_.reserve(2);

  energy_grid_.push_back(1.E-11);
  energy_grid_.push_back(200.);

  laws_.push_back(std::make_shared<Isotropic>());
  laws_.push_back(std::make_shared<Isotropic>());
}

AngleDistribution::AngleDistribution(const ACE& ace, int locb)
    : energy_grid_(), laws_() {
  // Locb must be >= 0! If locb == -1, it means that there is
  // no angular distribution for the reaction (must use product distribution)
  if (locb < 0) {
    std::string mssg = "Must have locb >= 0. Was provided with locb = " +
                       std::to_string(locb) + ".";
    throw PNDLException(mssg);
  }

  if (locb > 0) {
    // Set index
    std::size_t i = static_cast<std::size_t>(ace.AND() + locb - 1);

    // Get number of energies
    uint32_t NE = ace.xss<uint32_t>(i);

    // Read in energies
    energy_grid_ = ace.xss(i + 1, NE);

    if (!std::is_sorted(energy_grid_.begin(), energy_grid_.end())) {
      std::string mssg = "The energy grid must be sorted. Occurred at lobc = " +
                         std::to_string(locb) + ".";
      throw PNDLException(mssg);
    }

    // Get each table
    for (uint32_t j = 0; j < NE; j++) {
      int l = ace.xss<int>(i + 1 + NE + j);
      uint32_t loc = static_cast<uint32_t>(ace.AND() + std::abs(l) - 1);

      try {
        if (l > 0) {
          laws_.push_back(std::make_shared<EquiprobableAngleBins>(ace, loc));
        } else if (l < 0) {
          laws_.push_back(std::make_shared<AngleTable>(ace, loc));
        } else {
          laws_.push_back(std::make_shared<Isotropic>());
        }
      } catch (PNDLException& err) {
        std::string mssg =
            "Could not construct angular distribution for energy index = " +
            std::to_string(j) +
            ", energy = " + std::to_string(energy_grid_[j]) +
            " MeV. Occurred at loc = " + std::to_string(loc) +
            ", locb = " + std::to_string(locb) + ".";
        err.add_to_exception(mssg);
        throw err;
      }
    }
  } else if (locb == 0) {
    energy_grid_.push_back(1.E-5);
    laws_.push_back(std::make_shared<Isotropic>());
  }
}

AngleDistribution::AngleDistribution(
    const std::vector<double>& energy_grid,
    const std::vector<std::shared_ptr<AngleLaw>>& laws)
    : energy_grid_(energy_grid), laws_(laws) {
  if (energy_grid_.size() != laws_.size()) {
    std::string mssg =
        "The energy grid and the vector of laws must have the same size.";
    throw PNDLException(mssg);
  }

  if (!std::is_sorted(energy_grid_.begin(), energy_grid_.end())) {
    std::string mssg = "The energy grid must be sorted.";
    throw PNDLException(mssg);
  }
}

double AngleDistribution::sample_angle(
    double E_in, const std::function<double()>& rng) const {
  auto E_it = std::lower_bound(energy_grid_.begin(), energy_grid_.end(), E_in);
  if (E_it == energy_grid_.begin())
    return laws_.front()->sample_mu(rng);
  else if (E_it == energy_grid_.end())
    return laws_.back()->sample_mu(rng);
  E_it--;

  // Get index of low energy
  std::size_t l =
      static_cast<std::size_t>(std::distance(energy_grid_.begin(), E_it));
  double f = (E_in - energy_grid_[l]) / (energy_grid_[l + 1] - energy_grid_[l]);

  double mu = 0;

  if (rng() > f)
    mu = laws_[l]->sample_mu(rng);
  else
    mu = laws_[l + 1]->sample_mu(rng);

  if (std::abs(mu) > 1.) mu = std::copysign(1., mu);

  return mu;
}

double AngleDistribution::pdf(double E_in, double mu) const {
  auto E_it = std::lower_bound(energy_grid_.begin(), energy_grid_.end(), E_in);
  if (E_it == energy_grid_.begin())
    return laws_.front()->pdf(mu);
  else if (E_it == energy_grid_.end())
    return laws_.back()->pdf(mu);
  E_it--;

  // Get index of low energy
  std::size_t l =
      static_cast<std::size_t>(std::distance(energy_grid_.begin(), E_it));
  double f = (E_in - energy_grid_[l]) / (energy_grid_[l + 1] - energy_grid_[l]);

  return (1. - f) * laws_[l]->pdf(mu) + f * laws_[l + 1]->pdf(mu);
}

}  // namespace pndl
