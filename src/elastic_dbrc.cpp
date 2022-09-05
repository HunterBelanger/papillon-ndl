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

#include <PapillonNDL/elastic_dbrc.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <cmath>
#include <functional>

#include "constants.hpp"
#include "svt.hpp"
#include "vector.hpp"

namespace pndl {

ElasticDBRC::ElasticDBRC(const CrossSection& xs, const AngleDistribution& angle,
                         double awr, double temperature, bool use_tar,
                         double tar_threshold)
    : xs_(xs),
      angle_(angle),
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

double ElasticDBRC::temperature() const { return kT_ * MEV_TO_EV * EV_TO_K; }

double find_max_xs_value(const CrossSection& xs, const double& Emin,
                         const double& Emax) {
  const std::size_t i_min = xs.energy_grid().get_lower_index(Emin);
  // const double xs_Emin = xs(Emin, i_min);
  const double xs_Emin = xs(Emin);
  const std::size_t i_max = xs.energy_grid().get_lower_index(Emax);
  // const double xs_Emax = xs(Emax, i_max);
  const double xs_Emax = xs(Emax);
  double xs_max = std::max(xs_Emin, xs_Emax);

  for (std::size_t i = i_min + 1; i <= i_max; i++) {
    if (xs.xs()[i] > xs_max) xs_max = xs.xs()[i];
  }

  return xs_max;
}

Vector sample_target_velocity_dbrc(const double& Ein, const CrossSection& xs,
                                   const double& kT, const double& awr,
                                   const std::function<double()>& rng) {
  // Get min and max energies for finding the max, based on incident energy
  const Vector v_n(0., 0., std::sqrt(Ein));
  const double y = std::sqrt(awr * Ein / kT);
  const double y_min = std::max(0., y - 4.);
  const double Er_min = y_min * y_min * kT / awr;
  const double y_max = y + 4.;
  const double Er_max = y_max * y_max * kT / awr;
  const double xs_max = find_max_xs_value(xs, Er_min, Er_max);

  Vector v_t(0., 0., 0.);
  bool sample_velocity = true;
  while (sample_velocity) {
    v_t = sample_target_velocity(Ein, kT, awr, rng);

    const Vector vr = v_n - v_t;
    const double Er = vr.dot(vr);

    if (Er < Er_min || Er > Er_max) continue;

    const std::size_t i_Er = xs.energy_grid().get_lower_index(Er);
    const double xs_Er = xs(Er, i_Er);

    if (rng() * xs_max < xs_Er) {
      sample_velocity = false;
    }
  }

  return v_t;
}

AngleEnergyPacket ElasticDBRC::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
  // Direction in
  const Vector u_n(0., 0., 1.);

  // Get the "velocity" of the incident neutron in the lab frame
  const Vector v_n = u_n * std::sqrt(E_in);

  // Get the "velocity" of the target nuclide
  const bool TAR = (use_tar_ && E_in >= tar_threshold_ * kT_) ? true : false;
  const Vector v_t =
      TAR ? Vector(0., 0., 0.)
          : sample_target_velocity_dbrc(E_in, xs_, kT_, awr_, rng);

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
