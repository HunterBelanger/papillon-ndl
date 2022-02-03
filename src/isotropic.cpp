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
#include <PapillonNDL/isotropic.hpp>
#include <cmath>
#include <functional>

namespace pndl {

double Isotropic::sample_mu(std::function<double()> rng) const {
  double mu = 2. * rng() - 1.;
  if (std::abs(mu) > 1.) mu = std::copysign(1., mu);
  return mu;
}

double Isotropic::pdf(double mu) const {
  if (mu < -1. || mu > 1.) return 0.;
  return 0.5;
}

}  // namespace pndl
