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
#include <PapillonNDL/angle_table.hpp>
#include <PapillonNDL/pndl_exception.hpp>

namespace pndl {

AngleTable::AngleTable(const ACE& ace, std::size_t i) : distribution_(ace, i) {
  if (distribution_.min_value() < -1.) {
    std::string mssg =
        "Lowest posible cosine value is -1. Lowest given cosine is " +
        std::to_string(distribution_.min_value()) +
        ". Index to XSS block for table is " + std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  if (distribution_.max_value() > 1.) {
    std::string mssg =
        "Largest posible cosine value is 1. Largest given cosine is " +
        std::to_string(distribution_.max_value()) +
        ". Index to XSS block for table is " + std::to_string(i) + ".";
    throw PNDLException(mssg);
  }
}

AngleTable::AngleTable(const std::vector<double>& cosines,
                       const std::vector<double>& pdf,
                       const std::vector<double>& cdf, Interpolation interp)
    : distribution_(cosines, pdf, cdf, interp) {
  if (distribution_.min_value() < -1.) {
    std::string mssg =
        "Lowest posible cosine value is -1. Lowest given cosine is " +
        std::to_string(distribution_.min_value()) + ".";
    throw PNDLException(mssg);
  }

  if (distribution_.max_value() > 1.) {
    std::string mssg =
        "Largest posible cosine value is 1. Largest given cosine is " +
        std::to_string(distribution_.max_value()) + ".";
    throw PNDLException(mssg);
  }
}

AngleTable::AngleTable(const PCTable& table) : distribution_(table) {
  if (distribution_.min_value() < -1.) {
    std::string mssg =
        "Lowest posible cosine value is -1. Lowest given cosine is " +
        std::to_string(distribution_.min_value()) + ".";
    throw PNDLException(mssg);
  }

  if (distribution_.max_value() > 1.) {
    std::string mssg =
        "Largest posible cosine value is 1. Largest given cosine is " +
        std::to_string(distribution_.max_value()) + ".";
    throw PNDLException(mssg);
  }
}

double AngleTable::sample_mu(std::function<double()> rng) const {
  double mu = distribution_.sample_value(rng());
  if (std::abs(mu) > 1.) mu = std::copysign(1., mu);
  return mu;
}

}  // namespace pndl
