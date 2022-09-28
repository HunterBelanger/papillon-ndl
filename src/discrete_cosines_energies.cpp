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
#include <PapillonNDL/discrete_cosines_energies.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>
#include <optional>

namespace pndl {

// DiscreteCosinesEnergies can only be used with STIncoherentInelastic,
// and is the only distribution given, so the probability is always 1.
DiscreteCosinesEnergies::DiscreteCosinesEnergies(const ACE& ace)
    : incoming_energy_(), outgoing_energies_(), Noe(0), Nmu(0), skewed_(false) {
  // Make sure the distributions is discrete cosines and energies
  int32_t nxs_7 = ace.nxs(6);

  if (nxs_7 != 0 && nxs_7 != 1) {
    std::string mssg =
        "The provided ACE file does not contain a distribution of this form "
        "for Incoherent Inelastic scattering.";
    throw PNDLException(mssg);
  }

  if (nxs_7 == 1) {
    skewed_ = true;
  }

  // Read incident energy grid
  int32_t S = ace.jxs(0) - 1;
  uint32_t Ne = ace.xss<uint32_t>(S);  // Number of grid points
  incoming_energy_ = ace.xss(S + 1, Ne);

  // Make sure incident energy grid is sorted
  if (!std::is_sorted(incoming_energy_.begin(), incoming_energy_.end())) {
    std::string mssg = "The incident energy grid is not sorted.";
    throw PNDLException(mssg);
  }

  // Get the number of outgoing discrete energies
  Noe = static_cast<uint32_t>(ace.nxs(3));
  if (skewed_ && Noe < 6) {
    std::string mssg =
        "A skewed distribution must have at least 6 outgoing energies. Only " +
        std::to_string(Noe) + " were provided.";
    throw PNDLException(mssg);
  }

  // Get the number of outgoing discrete cosines
  Nmu = static_cast<uint32_t>(ace.nxs(2)) + 1;

  // Get the starting index for the distribution data
  int32_t i = ace.jxs(2) - 1;

  // Go through all incident energies
  for (std::size_t ie = 0; ie < Ne; ie++) {
    // Add the empy vector of all discrete energies
    outgoing_energies_.push_back({});

    // Get all outgoing energies for the current incident energy
    for (std::size_t oe = 0; oe < Noe; oe++) {
      double E_out = ace.xss(i);
      // Check E_out
      if (E_out <= 0.) {
        std::string mssg = "Nevative outgoing energy found at index " +
                           std::to_string(i) + ".";
        throw PNDLException(mssg);
      }
      i++;

      std::vector<double> mu = ace.xss(i, Nmu);
      // Check mu grid
      for (std::size_t j = 0; j < Nmu; j++) {
        if (mu[j] < -1. || mu[j] > 1.) {
          std::string mssg = "Invalid cosine value found at index " +
                             std::to_string(i + j) + ".";
          throw PNDLException(mssg);
        }
      }
      i += Nmu;

      // Add the pair to the last outgoing_energy
      outgoing_energies_.back().push_back({E_out, mu});
    }  // For all outgoing energies
  }    // For all incident energies
}

AngleEnergyPacket DiscreteCosinesEnergies::sample_angle_energy(
    double E_in, const std::function<double()>& rng) const {
  uint32_t j = 0;
  uint32_t k = static_cast<uint32_t>(Nmu * rng());
  // Sample j for the outgoing energy point
  if (!skewed_) {
    j = static_cast<uint32_t>(Noe * rng());
  } else {
    // Unit probability so that the sum of the probability of all bins
    // is 1 under the skewed conditions.
    double c = 1. / (10. * (static_cast<double>(Noe) - 3.));

    // Sample random value
    double xi = rng();

    if (xi >= 5 * c && xi < 1. - 5 * c) {
      // Select any one of the inner values, all with equal probability.
      // We put the most probable first to make less comparisons.
      j = static_cast<uint32_t>((Noe - 4) * rng() + 2);
    } else if (xi < c) {
      j = 0;
    } else if (xi >= c && xi < 5 * c) {
      j = 1;
    } else if (xi >= 1. - 5 * c && xi < 1. - c) {
      j = Noe - 2;
    } else {
      j = Noe - 1;
    }
  }

  // Find the incoming energy grid point
  auto Eit =
      std::lower_bound(incoming_energy_.begin(), incoming_energy_.end(), E_in);
  if (Eit == incoming_energy_.begin()) {
    // Below lowest incident energy.
    return {outgoing_energies_.front()[j].cosines[k],
            outgoing_energies_.front()[j].energy};
  } else if (Eit == incoming_energy_.end()) {
    // Above largest incident energy.
    return {outgoing_energies_.back()[j].cosines[k],
            outgoing_energies_.back()[j].energy};
  }
  std::size_t i = std::distance(incoming_energy_.begin(), Eit) - 1;
  double f = (E_in - incoming_energy_[i]) /
             (incoming_energy_[i + 1] - incoming_energy_[i]);

  double E_i_j = outgoing_energies_[i][j].energy;
  double E_i_1_j = outgoing_energies_[i + 1][j].energy;
  double E = E_i_j + f * (E_i_1_j - E_i_j);

  double mu_i_j_k = outgoing_energies_[i][j].cosines[k];
  double mu_i_1_j_k = outgoing_energies_[i + 1][j].cosines[k];
  double mu = mu_i_j_k + f * (mu_i_1_j_k - mu_i_j_k);

  if (std::abs(mu) > 1.) mu = std::copysign(1., mu);

  return {mu, E};
}

std::optional<double> DiscreteCosinesEnergies::angle_pdf(double /*E_in*/,
                                                         double /*mu*/) const {
  return std::nullopt;
}

std::optional<double> DiscreteCosinesEnergies::pdf(double /*E_in*/,
                                                   double /*mu*/,
                                                   double /*E_out*/) const {
  return std::nullopt;
}

}  // namespace pndl
