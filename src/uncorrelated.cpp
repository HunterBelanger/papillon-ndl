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
#include <PapillonNDL/uncorrelated.hpp>

namespace pndl {

Uncorrelated::Uncorrelated(const AngleDistribution& angle,
                           std::shared_ptr<EnergyLaw> energy)
    : angle_(angle), energy_(energy) {
  if (!energy_) {
    std::string mssg = "Provided energy distribution is nullptr.";
    throw PNDLException(mssg);
  }
}

AngleEnergyPacket Uncorrelated::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
  double mu = angle_.sample_angle(E_in, rng);
  double E_out = energy_->sample_energy(E_in, rng);
  return {mu, E_out};
}

std::optional<double> Uncorrelated::angle_pdf(double E_in, double mu) const {
  return angle_.pdf(E_in, mu);
}

std::optional<double> Uncorrelated::pdf(double E_in, double mu,
                                        double E_out) const {
  std::optional<double> energy_pdf = energy_->pdf(E_in, E_out);

  if (!energy_pdf) {
    return std::nullopt;
  }

  double angle_pdf = angle_.pdf(E_in, mu);

  energy_pdf.value() *= angle_pdf;

  return energy_pdf;
}

}  // namespace pndl
