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

/**
 * @file
 * @author Hunter Belanger
 */

#include "incoherent_elastic.hpp"
#include "constants.hpp"
#include "interpolator.hpp"

#include <Log.hpp>
using namespace njoy;

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <variant>

IncoherentElastic::IncoherentElastic(
    const section::Type<7, 2>::IncoherentElastic& ie)
    : W_(makeTab1(ie.boundaries(), ie.interpolants(), ie.temperatures(),
                  ie.debyeWallerValues())),
      bound_xs_(ie.boundCrossSection()) {

  // Write Information
  Log::info("Incoherent Elastic");
  Log::info("------------------");
  Log::info("Bound XS = {}", bound_xs_);
  Log::info("");
}

double IncoherentElastic::dxs(double T, double Ein, double mu) const {
  if (mu < -1. || mu > 1.) {
    throw std::runtime_error(
        "IncoherentElastic::dxs: mu must be in interval [-1, 1].");
  }

  return (bound_xs_ / (4. * PI)) * std::exp(-2. * Ein * W_(T) * (1. - mu));
}

double IncoherentElastic::xs(double T, double Ein) const {
  const double W = W_(T);
  return 0.5 * bound_xs_ * (1. - std::exp(-4. * Ein * W)) / (2. * Ein * W);
}
