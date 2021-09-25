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
#include <PapillonNDL/multiple_distribution.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <optional>

namespace pndl {

MultipleDistribution::MultipleDistribution(
    const std::vector<std::shared_ptr<AngleEnergy>>& distributions,
    const std::vector<std::shared_ptr<Tabulated1D>>& probabilities)
    : distributions_(distributions), probabilities_(probabilities) {
  if (distributions_.size() != probabilities_.size()) {
    std::string mssg =
        "Different number of distributions and probabilities provided.";
    throw PNDLException(mssg);
  }

  if (distributions_.size() <= 1) {
    std::string mssg = "At least two distributions must be provided.";
    throw PNDLException(mssg);
  }

  for (std::size_t i = 0; i < distributions_.size(); i++) {
    if (!distributions_[i]) {
      std::string mssg =
          "Distribution at index " + std::to_string(i) + " is nullptr.";
      throw PNDLException(mssg);
    } else if (!probabilities_[i]) {
      std::string mssg =
          "Probability at index " + std::to_string(i) + " is nullptr.";
      throw PNDLException(mssg);
    }
  }
}

AngleEnergyPacket MultipleDistribution::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
  // First select distribution
  double xi = rng();
  double sum = 0.;
  for (std::size_t d = 0; d < distributions_.size(); d++) {
    sum += (*probabilities_[d])(E_in);
    if (xi < sum) {
      return distributions_[d]->sample_angle_energy(E_in, rng);
    }
  }

  // Shouldn't get here, but if we do, use the last distribution
  return distributions_.back()->sample_angle_energy(E_in, rng);
}

std::optional<double> MultipleDistribution::angle_pdf(double E_in,
                                                      double mu) const {
  double a_pdf = 0.;
  for (std::size_t d = 0; d < distributions_.size(); d++) {
    std::optional<double> E_in_pdf = distributions_[d]->angle_pdf(E_in, mu);

    if (!E_in_pdf) {
      return std::nullopt;
    }

    a_pdf += (*probabilities_[d])(E_in)*E_in_pdf.value();
  }

  return a_pdf;
}

std::optional<double> MultipleDistribution::pdf(double E_in, double mu,
                                                double E_out) const {
  double j_pdf = 0.;
  for (std::size_t d = 0; d < distributions_.size(); d++) {
    std::optional<double> E_in_pdf = distributions_[d]->pdf(E_in, mu, E_out);

    if (!E_in_pdf) {
      return std::nullopt;
    }

    j_pdf += (*probabilities_[d])(E_in)*E_in_pdf.value();
  }

  return j_pdf;
}

}  // namespace pndl
