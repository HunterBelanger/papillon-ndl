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

#include <PapillonNDL/elastic_dbrc.hpp>
#include <PapillonNDL/elastic_svt.hpp>
#include <cmath>
#include <functional>

#include "vector.hpp"

namespace pndl {

double ElasticDBRC::max_xs_value(const double& Emin, const double& Emax) const {
  const std::size_t i_min = xs_.energy_grid().get_lower_index(Emin);
  // const double xs_Emin = xs(Emin, i_min);
  const double xs_Emin = xs_(Emin);
  const std::size_t i_max = xs_.energy_grid().get_lower_index(Emax);
  // const double xs_Emax = xs(Emax, i_max);
  const double xs_Emax = xs_(Emax);
  double xs_max = std::max(xs_Emin, xs_Emax);

  for (std::size_t i = i_min + 1; i <= i_max; i++) {
    if (xs_.xs()[i] > xs_max) xs_max = xs_.xs()[i];
  }

  return xs_max;
}

std::array<double, 3> ElasticDBRC::sample_target_velocity(
    const double& Ein, const double& kT, const double& awr,
    const std::function<double()>& rng) const {
  const static ElasticSVT svt;

  // Get min and max energies for finding the max, based on incident energy
  const Vector v_n(0., 0., std::sqrt(Ein));
  const double y = std::sqrt(awr * Ein / kT);
  const double y_min = std::max(0., y - 4.);
  const double Er_min = y_min * y_min * kT / awr;
  const double y_max = y + 4.;
  const double Er_max = y_max * y_max * kT / awr;
  const double xs_max = this->max_xs_value(Er_min, Er_max);

  Vector v_t(0., 0., 0.);
  bool sample_velocity = true;
  while (sample_velocity) {
    v_t = svt.sample_target_velocity(Ein, kT, awr, rng);

    const Vector vr = v_n - v_t;
    const double Er = vr.dot(vr);

    if (Er < Er_min || Er > Er_max) continue;

    const std::size_t i_Er = xs_.energy_grid().get_lower_index(Er);
    const double xs_Er = xs_(Er, i_Er);

    if (rng() * xs_max < xs_Er) {
      sample_velocity = false;
    }
  }

  return v_t.array();
}

std::string ElasticDBRC::algorithm() const { return "DBRC"; }

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
