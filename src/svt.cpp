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

}  // namespace pndl
