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
#include <PapillonNDL/cm_distribution.hpp>
#include <PapillonNDL/frame.hpp>

/**
 * @file
 * @author Hunter Belanger
 */

namespace pndl {

AngleEnergyPacket CMDistribution::sample_angle_energy(
    double E_in, const std::function<double()>& rng) const {
  AngleEnergyPacket out = distribution_->sample_angle_energy(E_in, rng);

  CMToLab::transform(E_in, awr_, out);

  return out;
}

std::optional<double> CMDistribution::angle_pdf(double E_in, double mu) const {
  // First we need the angle in the CM frame
  auto cm_angles = LabToCM::angle(E_in, awr_, q_, mu);

  // There can be at most up to two angles in the CM frame for a given
  // angle in the Lab frame. If one angle is returned, then we only
  // return the PDF component for that angle. If two angles are returned,
  // we take the sum of their respective PDFs. If no angles are returned,
  // this means that it is imposible to have a scattering angle of mu in
  // the lab frame for the given reaction. In this case, zero is returned.

  // Variables to contain the components
  double p1 = 0.;
  double p2 = 0.;

  // Treat first angle
  if (cm_angles.first) {
    double mu_cm = cm_angles.first.value();
    auto p1_opt = distribution_->angle_pdf(E_in, mu_cm);
    if (p1_opt) {
      p1 = p1_opt.value() * CMToLab::angle_jacobian(E_in, awr_, q_, mu, mu_cm);
    }
  }

  // Treat second angle
  if (cm_angles.second) {
    double mu_cm = cm_angles.second.value();
    auto p2_opt = distribution_->angle_pdf(E_in, mu_cm);
    if (p2_opt) {
      p2 = p2_opt.value() * CMToLab::angle_jacobian(E_in, awr_, q_, mu, mu_cm);
    }
  }

  return p1 + p2;
}

std::optional<double> CMDistribution::pdf(double E_in, double mu,
                                          double E_out) const {
  // First, we need to get the angle and energy in the CM frame
  double mu_cm = mu;
  double Eout_cm = E_out;
  LabToCM::transform(E_in, awr_, mu_cm, Eout_cm);

  // We now get the PDF in the CM frame
  auto p = distribution_->pdf(E_in, mu_cm, Eout_cm);

  // Convert the CM frame to the lab frame
  if (p) {
    p.value() *= CMToLab::jacobian(E_out, Eout_cm);
  }

  return p;
}

}  // namespace pndl
