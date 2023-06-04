/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2023, Hunter Belanger
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

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/direct_sab.hpp>
#include <PapillonNDL/interpolation.hpp>
#include <PapillonNDL/pctable.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <cmath>
#include <sstream>
#include <vector>

#include "constants.hpp"

namespace pndl {

DirectSab::DirectSab(const ACE& ace)
    : AngleEnergy(), incoming_energy_(), beta_dists_(), kT_(0.), A_(0.) {
  // Make sure this is a special TSL ACE made by Panglos
  const int32_t nxs_7 = ace.nxs(6);
  if (nxs_7 != 3) {
    std::string mssg =
        "ACE File does not contain a Direct S(a,b) distribution.";
    throw PNDLException(mssg);
  }

  // Read the incident energy grid
  try {
    const std::size_t S = static_cast<std::size_t>(ace.jxs(0) - 1);
    const uint32_t Ne = ace.xss<uint32_t>(S);  // Number of grid points
    incoming_energy_ = ace.xss(S + 1, Ne);
  } catch (PNDLException& err) {
    std::string mssg = "Could not read incoming energy grid.";
    err.add_to_exception(mssg);
    throw err;
  }

  const std::size_t BPS = static_cast<std::size_t>(
      ace.jxs(2) - 1);  // Starting location of the beta pointers
  beta_dists_.reserve(incoming_energy_.size());
  for (std::size_t i = 0; i < incoming_energy_.size(); i++) {
    const std::size_t BLOC = static_cast<std::size_t>(ace.jxs(1)) +
                             ace.xss<std::size_t>(BPS + i) - 1;
    const std::size_t NB = ace.xss<std::size_t>(BLOC);
    const std::vector<double> beta = ace.xss(BLOC + 1, NB);
    const std::vector<double> pdf = ace.xss(BLOC + 1 + NB, NB);
    const std::vector<double> cdf = ace.xss(BLOC + 1 + NB + NB, NB);

    std::vector<PCTable> alpha_dists;
    alpha_dists.reserve(beta.size());
    const std::size_t APS = BLOC + 1 + NB + NB + NB;
    for (std::size_t b = 0; b < beta.size(); b++) {
      const std::size_t ALOC = static_cast<std::size_t>(ace.jxs(1)) +
                               ace.xss<std::size_t>(APS + b) - 1;
      ;
      const std::size_t NA = ace.xss<std::size_t>(ALOC);
      const std::vector<double> alpha = ace.xss(ALOC + 1, NA);
      const std::vector<double> apdf = ace.xss(ALOC + 1 + NA, NA);
      const std::vector<double> acdf = ace.xss(ALOC + 1 + NA + NA, NA);
      try {
        alpha_dists.emplace_back(alpha, apdf, acdf, Interpolation::LinLin);
      } catch (PNDLException& err) {
        std::stringstream mssg;
        mssg << "Could not create alpha distribution for incident energy index"
             << i << ", beta index " << b << ".";
        err.add_to_exception(mssg.str());
        throw err;
      }
    }

    try {
      beta_dists_.emplace_back(beta, pdf, cdf, alpha_dists);
    } catch (PNDLException& err) {
      std::stringstream mssg;
      mssg << "Could not create BetaAlphaTable for incident energy index " << i
           << ".";
      err.add_to_exception(mssg.str());
      throw err;
    }
  }

  // Get the AWR
  A_ = ace.awr();

  // Get the temperature
  kT_ = ace.temperature() * K_TO_EV * EV_TO_MEV;
}

AngleEnergyPacket DirectSab::sample_angle_energy(
    double E_in, const std::function<double()>& rng) const {
  // Determine the index of the bounding tabulated incoming energies
  std::size_t l;
  double f;  // Interpolation factor
  auto in_E_it =
      std::lower_bound(incoming_energy_.begin(), incoming_energy_.end(), E_in);
  if (in_E_it == incoming_energy_.begin()) {
    l = 0;
    f = 0.;
  } else if (in_E_it == incoming_energy_.end()) {
    l = incoming_energy_.size() - 2;
    f = 1.;
  } else {
    l = static_cast<std::size_t>(
        std::distance(incoming_energy_.begin(), in_E_it) - 1);
    f = (E_in - incoming_energy_[l]) /
        (incoming_energy_[l + 1] - incoming_energy_[l]);
  }

  // Sample outgoing energy, and get R and A for mu
  const double B_i_M = beta_dists_[l].max_beta();
  const double B_i_1_M = beta_dists_[l + 1].max_beta();
  const double Bmin = -E_in / kT_;
  const double Bmax = B_i_M + f * (B_i_1_M - B_i_M);

  double E_out = 0.;
  double mu = 0.;
  bool sampled = false;
  while (sampled == false) {
    AlphaBetaPacket tmp{0., 0.};
    double B_l_1{0.}, B_l_M{0.};
    if (rng() > f) {
      tmp = beta_dists_[l].sample_alpha_beta(rng);
      B_l_1 = beta_dists_[l].min_beta();
      B_l_M = beta_dists_[l].max_beta();
    } else {
      tmp = beta_dists_[l + 1].sample_alpha_beta(rng);
      B_l_1 = beta_dists_[l + 1].min_beta();
      B_l_M = beta_dists_[l + 1].max_beta();
    }

    // Cacluate outgoing energy
    const double B =
        Bmin + ((tmp.beta - B_l_1) / (B_l_M - B_l_1)) * (Bmax - Bmin);
    E_out = B * kT_ + E_in;
    if (E_out < 0.) E_out = 0.;

    // Calculate scattering angle
    mu = (E_in + E_out - tmp.alpha * A_ * kT_) / (2. * std::sqrt(E_in * E_out));
    if (std::abs(mu) <= 1.) sampled = true;
  }

  return {mu, E_out};
}

double DirectSab::temperature() const { return kT_ * MEV_TO_EV * EV_TO_K; }

}  // namespace pndl