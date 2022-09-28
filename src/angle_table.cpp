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
#include <PapillonNDL/angle_table.hpp>
#include <PapillonNDL/linearize.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <functional>

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

AngleTable::AngleTable(const Legendre& legendre)
    : distribution_({-1., 1.}, {0.5, 0.5}, {0., 1.}, Interpolation::LinLin) {
  std::function<double(double)> l =
      std::bind(&Legendre::pdf, legendre, std::placeholders::_1);
  try {
    Tabulated1D tab_pdf = linearize(-1., 1., l);

    std::vector<double> cosines = tab_pdf.x();
    std::vector<double> pdf = tab_pdf.y();
    std::vector<double> cdf(pdf.size(), 0.);

    // Trapezoid rule
    for (std::size_t i = 0; i < cosines.size() - 1; i++) {
      cdf[i + 1] =
          0.5 * (pdf[i] + pdf[i + 1]) * (cosines[i + 1] - cosines[i]) + cdf[i];
    }

    const double norm = cdf.back();
    for (std::size_t i = 0; i < cosines.size(); i++) {
      pdf[i] /= norm;
      cdf[i] /= norm;
    }

    distribution_ = PCTable(cosines, pdf, cdf, Interpolation::LinLin);
  } catch (PNDLException& err) {
    std::string mssg = "Could not linearize Legenre distribution.";
    err.add_to_exception(mssg);
    throw err;
  }
}

double AngleTable::sample_mu(const std::function<double()>& rng) const {
  double mu = distribution_.sample_value(rng());
  if (std::abs(mu) > 1.) mu = std::copysign(1., mu);
  return mu;
}

}  // namespace pndl
