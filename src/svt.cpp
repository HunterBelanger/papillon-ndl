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

#include "svt.hpp"

#include <cmath>

namespace pndl {

Vector sample_target_velocity(const double& Ein, const double& kT,
                              const double& awr,
                              const std::function<double()>& rng) {
  const double y = std::sqrt(awr * Ein / kT);
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
  double s_t = std::sqrt(x_sqrd * kT / awr);

  // Use mu to get the direction vector of the target. We know in the sample
  // method that we always assume the same incident neutron vector:
  const Vector u_n{0., 0., 1.};
  const Vector u_t = u_n.rotate(mu, 2. * PI * rng());

  return u_t * s_t;
}

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

}  // namespace pndl
