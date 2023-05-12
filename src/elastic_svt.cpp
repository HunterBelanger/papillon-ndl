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

#include <PapillonNDL/elastic_svt.hpp>
#include <cmath>
#include <functional>

#include "vector.hpp"

namespace pndl {

std::array<double, 3> ElasticSVT::sample_target_velocity(
    const double& Ein, const double& kT, const double& awr,
    const std::function<double()>& rng) const {
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

  return (u_t * s_t).array();
}

std::string ElasticSVT::algorithm() const { return "SVT"; }

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
