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

#include <PapillonNDL/elastic.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <cmath>

#include "constants.hpp"
#include "vector.hpp"

namespace pndl {

Elastic::Elastic(std::shared_ptr<ElasticDopplerBroadener> broadener,
                 const AngleDistribution& angle, double awr, double temperature,
                 bool use_tar, double tar_threshold)
    : angle_(angle),
      awr_(awr),
      kT_(temperature * K_TO_EV * EV_TO_MEV),
      use_tar_(use_tar),
      tar_threshold_(tar_threshold),
      broadener_(broadener) {
  if (awr_ <= 0.) {
    std::string mssg = "Atomic weight ratio must be greater than zero.";
    throw PNDLException(mssg);
  }

  if (kT_ < 0.) {
    std::string mssg = "Temperature must be greater than or equal to zero.";
    throw PNDLException(mssg);
  }

  if (tar_threshold_ < 0.) {
    std::string mssg =
        "Target At Rest threshold must be greater than or equal to zero.";
    throw PNDLException(mssg);
  }

  if (broadener_ == nullptr) {
    std::string mssg = "Provided ElasticDopplerBroadener pointer was nullptr.";
    throw PNDLException(mssg);
  }
}

AngleEnergyPacket Elastic::sample_angle_energy(
    double E_in, const std::function<double()>& rng) const {
  // Direction in
  const Vector u_n(0., 0., 1.);

  // Get the "velocity" of the incident neutron in the lab frame
  const Vector v_n = u_n * std::sqrt(E_in);

  // Get the "velocity" of the target nuclide
  const bool TAR = (use_tar_ && E_in >= tar_threshold_ * kT_) ? true : false;
  const Vector v_t =
      TAR ? Vector(0., 0., 0.)
          : broadener_->sample_target_velocity(E_in, kT_, awr_, rng);

  // Calculate the "velocity" of the center of mass.
  const Vector v_cm = (v_n + (v_t * awr_)) / (awr_ + 1.);

  // Calculate the "velocity" of the neutron in the CM frame
  const Vector V_n = v_n - v_cm;

  // Calculate the "speed" in the CM frame
  const double S_n = V_n.magnitude();

  // Calculate direction in the CM frame
  const Vector U_n = V_n / S_n;

  // Sample the scattering angle in the CM frame
  const double mu_cm = angle_.sample_angle(E_in, rng);

  // Get the outgoing velocity in the CM frame
  const Vector V_n_out = U_n.rotate(mu_cm, 2. * PI * rng()) * S_n;

  // Get the outgoing velocity in the LAB frame
  const Vector v_n_out = V_n_out + v_cm;

  // Calculate "speed" in the LAB frame
  const double s_n_out = v_n_out.magnitude();

  // Calculate direction in the LAB frrame
  const Vector u_n_out = v_n_out / s_n_out;

  // Calculate outgoing energy in LAB frame
  const double E_out = s_n_out * s_n_out;

  // Calculate the cosine of the scattering angle in the LAB frame
  const double mu_lab = u_n.dot(u_n_out);

  return {mu_lab, E_out};
}

void Elastic::set_elastic_doppler_broadener(
    std::shared_ptr<ElasticDopplerBroadener> broadener) {
  if (broadener == nullptr) {
    std::string mssg = "Provided ElasticDopplerBroadener pointer was nullptr.";
    throw PNDLException(mssg);
  }

  broadener_ = broadener;
}

double Elastic::temperature() const { return kT_ * MEV_TO_EV * EV_TO_K; }

void Elastic::set_temperature(double temperature) {
  kT_ = temperature * K_TO_EV * EV_TO_MEV;
}

void Elastic::set_tar_threshold(double tar_threshold) {
  tar_threshold_ = tar_threshold;

  if (tar_threshold_ < 0.) {
    std::string mssg =
        "Target At Rest threshold must be greater than or equal to zero.";
    throw PNDLException(mssg);
  }
}

}  // namespace pndl
