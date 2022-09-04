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
#ifndef PAPILLON_NDL_VECTOR_H
#define PAPILLON_NDL_VECTOR_H

#include <cmath>

#include "constants.hpp"

namespace pndl {
struct Vector {
  double x, y, z;

  Vector(double x, double y, double z) : x(x), y(y), z(z) {}

  Vector operator+(const Vector& v) const {
    return {x + v.x, y + v.y, z + v.z};
  }

  Vector operator-(const Vector& v) const {
    return {x - v.x, y - v.y, z - v.z};
  }

  Vector operator*(const double& c) const { return {x * c, y * c, z * c}; }

  Vector operator/(const double& c) const { return {x / c, y / c, z / c}; }

  double dot(const Vector& v) const { return x * v.x + y * v.y + z * v.z; }

  double magnitude() const { return std::sqrt(x * x + y * y + z * z); }

  Vector rotate(double mu, double phi) const {
    if (mu < -1.)
      mu = -1.;
    else if (mu > 1.)
      mu = 1.;

    if (phi < 0.)
      phi = 0.;
    else if (phi > 2. * PI)
      phi = 2. * PI;

    double xo, yo, zo;
    const double c = std::cos(phi);
    const double s = std::sin(phi);
    const double C = std::sqrt(1. - mu * mu);

    if (std::abs(1. - z * z) > 1.E-10) {
      const double denom = std::sqrt(1. - z * z);

      xo = x * mu + C * (c * x * z - s * y) / denom;
      yo = y * mu + C * (c * y * z + s * x) / denom;
      zo = z * mu - c * C * denom;
    } else {
      const double denom = std::sqrt(1. - y * y);

      xo = x * mu + C * (c * x * y + s * z) / denom;
      yo = y * mu - c * C * denom;
      zo = z * mu + C * (c * y * z - s * x) / denom;
    }

    return {xo, yo, zo};
  }
};
}  // namespace pndl

#endif
