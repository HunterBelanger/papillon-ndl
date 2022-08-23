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

#include <PapillonNDL/elastic_svt.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <cmath>
#include <functional>

#include "constants.hpp"

namespace pndl {

ElasticSVT::ElasticSVT(const AngleDistribution& angle, double awr,
                       double temperature, double tar_threshold)
    : angle_(angle),
      awr_(awr),
      temperature_(temperature),
      tar_threshold_(tar_threshold) {
  if (awr_ <= 0.) {
    std::string mssg = "Atomic weight ratio must be greater than zero.";
    throw PNDLException(mssg);
  }

  if (temperature_ < 0.) {
    std::string mssg = "Temperature must be greater than or equal to zero.";
    throw PNDLException(mssg);
  }

  if (tar_threshold_ < 0.) {
    std::string mssg =
        "Target At Rest threshold must be greater than or equal to zero.";
    throw PNDLException(mssg);
  }
}

inline double ElasticSVT::Vector::magnitude() const {
  return std::sqrt(x * x + y * y + z * z);
}

inline ElasticSVT::Vector ElasticSVT::Vector::rotate(double mu, double phi) const {
  double xo, yo, zo;
  const double c = std::cos(phi);
  const double s = std::sin(phi);
  const double C = std::sqrt(1. - mu*mu);

  if (std::abs(1. - z*z) > 1.E-10) {
    const double denom = std::sqrt(1. - z*z);

    xo = x*mu + C*(c*x*z - s*y)/denom;
    yo = y*mu + C*(c*y*z + s*x)/denom;
    zo = z*mu - c*C*denom;
  } else {
    const double denom = std::sqrt(1. - y*y);

    xo = x*mu + C*(c*x*y + s*z)/denom;
    yo = y*mu - c*C*denom;
    zo = z*mu + C*(c*y*z - s*x)/denom;
  }

  return {xo, yo, zo};
}

AngleEnergyPacket ElasticSVT::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
  // Get the "velocity" of the incident neutron in the lab frame
  const Vector v_n = Vector(0., 0., 1.) * std::sqrt(E_in);

  // Get the "velocity" of the target nuclide
  const Vector v_t = sample_target_velocity(rng);

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
  const double S_out = v_n_out.magnitude();

  // Calculate outgoing energy in LAB frame
  const double E_out = S_out * S_out;

  // Calculate the cosine of the scattering angle in the LAB frame
  const double mu_lab = v_n.dot(v_n_out);

  return {mu_lab, E_out};
}

ElasticSVT::Vector ElasticSVT::sample_target_velocity(
    double Ein, std::function<double()> rng) const {
  if (Ein >= tar_threshold_ * temperature_ * K_TO_EV * EV_TO_MEV && awr_ > 1.) {
    return {0., 0., 0.};
  }

  const double y =
      std::sqrt(0.5 * awr_ * Ein / (temperature_ * K_TO_EV * EV_TO_MEV));
  double x_sqrd = 0.;
  double mu = 0.;

  const double P_C49 = 2. / (std::sqrt(PI) * y + 2.);

  bool sample_velocity = true;
  while (sample_velocity) {
    if (rng() < P_C49) {
      // Sample x from the distribution C49 in MC sampler
      x_sqrd = -std::log(rng() * rng());
    } else {
      // Sample x from the distribution C61 in MC sampler
      const double c = std::cos(PI / 2.0 * rng());
      x_sqrd = -std::log(rng()) - std::log(rng()) * c * c;
    }

    const double x = std::sqrt(x_sqrd);
    mu = 2. * rng() - 1.;
    const double P_accept =
        std::sqrt(y * y + x_sqrd - 2. * y * x * mu) / (x + y);

    if (rng() < P_accept) sample_velocity = false;
  }

  // Get speed of target
  double s_t =
      std::sqrt(x_sqrd * 2. * temperature_ * K_TO_EV * EV_TO_MEV / awr_);

  // Use mu to get the direction vector of the target. We know in the sample
  // method that we always assume the same incident neutron vector:
  const Vector u_n{0., 0., 1.};
  Vector u_t = u_n.rotate(mu, 2.*PI*rng());

  return u_t * s_t;
}

}  // namespace pndl
