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
#include "svt.hpp"
#include "vector.hpp"

namespace pndl {

ElasticSVT::ElasticSVT(const AngleDistribution& angle, double awr,
                       double temperature, bool use_tar, double tar_threshold)
    : angle_(angle),
      awr_(awr),
      kT_(temperature * K_TO_EV * EV_TO_MEV),
      use_tar_(use_tar),
      tar_threshold_(tar_threshold) {
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

  if (use_tar_ == false) tar_threshold_ = INF;
}

double ElasticSVT::temperature() const { return kT_ * MEV_TO_EV * EV_TO_K; }

AngleEnergyPacket ElasticSVT::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
  // Direction in
  const Vector u_n(0., 0., 1.);

  // Get the "velocity" of the incident neutron in the lab frame
  const Vector v_n = u_n * std::sqrt(E_in);

  // Get the "velocity" of the target nuclide
  const bool TAR = (use_tar_ && E_in >= tar_threshold_ * kT_) ? true : false;
  const Vector v_t =
      TAR ? Vector(0., 0., 0.) : sample_target_velocity(E_in, kT_, awr_, rng);

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

}  // namespace pndl

/*
 * REFERENCES
 *
 * [1] R. R. Coveyou, R. R. Bate, and R. K. Osborn, “Effect of moderator
 * temperature upon neutron flux in infinite, capturing medium,” J Nucl Energy
 * 1954, vol. 2, no. 3–4, pp. 153–167, 1956, doi: 10.1016/0891-3919(55)90030-9.
 *
 * [2] P. K. Romano and J. A. Walsh, “An improved target velocity sampling
 * algorithm for free gas elastic scattering,” Ann Nucl Energy, vol. 114, no.
 * Ann.  Nucl. Energy 36 2009, pp. 318–324, 2018,
 * doi: 10.1016/j.anucene.2017.12.044.
 */
