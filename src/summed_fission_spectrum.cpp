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

#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/summed_fission_spectrum.hpp>
#include <optional>

namespace pndl {

SummedFissionSpectrum::SummedFissionSpectrum(std::shared_ptr<STReaction> mt19,
                                             std::shared_ptr<STReaction> mt20,
                                             std::shared_ptr<STReaction> mt21,
                                             std::shared_ptr<STReaction> mt38)
    : reactions_{mt19, mt20, mt21, mt38} {
  if (mt19 == nullptr) {
    std::string mssg = "Reaction for MT=19 was nullptr.";
    throw PNDLException(mssg);
  } else if (mt20 == nullptr) {
    std::string mssg = "Reaction for MT=20 was nullptr.";
    throw PNDLException(mssg);
  } else if (mt21 == nullptr) {
    std::string mssg = "Reaction for MT=21 was nullptr.";
    throw PNDLException(mssg);
  } else if (mt38 == nullptr) {
    std::string mssg = "Reaction for MT=38 was nullptr.";
    throw PNDLException(mssg);
  }
}

inline void SummedFissionSpectrum::compute_probabilities(
    std::array<double, 4>& probs, double E_in) const {
  // Initialize arrays for the xs values and probabilities
  std::array<double, 4> xs{0., 0., 0., 0.};
  probs = {0., 0., 0., 0.};

  // Get the energy index from the energy grid
  std::size_t i = reactions_[0]->xs().energy_grid().get_lower_index(E_in);
  xs[0] = reactions_[0]->xs()(E_in, i);
  xs[1] = reactions_[1]->xs()(E_in, i);
  xs[2] = reactions_[2]->xs()(E_in, i);
  xs[3] = reactions_[3]->xs()(E_in, i);
  const double xs_tot = xs[0] + xs[1] + xs[2] + xs[3];

  // Calculate probabilities
  probs[0] = xs[0] / xs_tot;
  probs[1] = probs[0] + (xs[1] / xs_tot);
  probs[2] = probs[1] + (xs[2] / xs_tot);
  probs[3] = probs[2] + (xs[3] / xs_tot);
}

AngleEnergyPacket SummedFissionSpectrum::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
  // Get probabilities
  std::array<double, 4> probs{0., 0., 0., 0.};
  this->compute_probabilities(probs, E_in);

  // Sample the distribution to use
  const double xi = rng();
  if (xi < probs[0]) {
    return reactions_[0]->neutron_distribution().sample_angle_energy(E_in, rng);
  } else if (xi < probs[1]) {
    return reactions_[1]->neutron_distribution().sample_angle_energy(E_in, rng);
  } else if (xi < probs[2]) {
    return reactions_[2]->neutron_distribution().sample_angle_energy(E_in, rng);
  } else {
    return reactions_[3]->neutron_distribution().sample_angle_energy(E_in, rng);
  }
}

std::optional<double> SummedFissionSpectrum::angle_pdf(double E_in, double mu) const {
  // Get probabilities
  std::array<double, 4> probs {0., 0., 0., 0.};
  this->compute_probabilities(probs, E_in);

  double pdf = 0;

  if (probs[0] > 0.) {
    const auto pdf_0 = reactions_[0]->neutron_distribution().angle_pdf(E_in, mu);
    if (pdf_0) pdf += probs[0] * pdf_0.value();
    else return std::nullopt;
  } else if (probs[1] > 0.) {
    const auto pdf_1 = reactions_[1]->neutron_distribution().angle_pdf(E_in, mu);
    if (pdf_1) pdf += (probs[1] - probs[0]) * pdf_1.value();
    else return std::nullopt;
  } else if (probs[2] > 0.) {
    const auto pdf_2 = reactions_[2]->neutron_distribution().angle_pdf(E_in, mu);
    if (pdf_2) pdf += (probs[2] - probs[1]) * pdf_2.value();
    else return std::nullopt;
  } else if (probs[3] > 0.) {
    const auto pdf_3 = reactions_[3]->neutron_distribution().angle_pdf(E_in, mu);
    if (pdf_3) pdf += (probs[3] - probs[2]) * pdf_3.value();
    else return std::nullopt;
  }

  return pdf;
}

std::optional<double> SummedFissionSpectrum::pdf(double E_in, double mu,
                                                 double E_out) const {
  // Get probabilities
  std::array<double, 4> probs {0., 0., 0., 0.};
  this->compute_probabilities(probs, E_in);

  double pdf = 0;

  if (probs[0] > 0.) {
    const auto pdf_0 = reactions_[0]->neutron_distribution().pdf(E_in, mu, E_out);
    if (pdf_0) pdf += probs[0] * pdf_0.value();
    else return std::nullopt;
  } else if (probs[1] > 0.) {
    const auto pdf_1 = reactions_[1]->neutron_distribution().pdf(E_in, mu, E_out);
    if (pdf_1) pdf += (probs[1] - probs[0]) * pdf_1.value();
    else return std::nullopt;
  } else if (probs[2] > 0.) {
    const auto pdf_2 = reactions_[2]->neutron_distribution().pdf(E_in, mu, E_out);
    if (pdf_2) pdf += (probs[2] - probs[1]) * pdf_2.value();
    else return std::nullopt;
  } else if (probs[3] > 0.) {
    const auto pdf_3 = reactions_[3]->neutron_distribution().pdf(E_in, mu, E_out);
    if (pdf_3) pdf += (probs[3] - probs[2]) * pdf_3.value();
    else return std::nullopt;
  }

  return pdf;
}

}  // namespace pndl
