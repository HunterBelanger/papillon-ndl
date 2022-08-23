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
#ifndef PAPILLON_NDL_ELASTIC_SVT_H
#define PAPILLON_NDL_ELASTIC_SVT_H

#include <PapillonNDL/angle_distribution.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <functional>
#include <optional>

/**
 * @file
 * @author Hunter Belanger
 */

namespace pndl {

class ElasticSVT : public AngleEnergy {
 public:
  ElasticSVT(const AngleDistribution& angle, double awr, double temperature,
             double tar_threshold);

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  std::optional<double> angle_pdf(double E_in, double mu) const override final {
    return std::nullopt;
  }

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override final {
    return std::nullopt;
  }

  /**
   * @brief Returns the AngleDistribution which describes the distribution for
   *        the cosine of the scattering angle in the center-of-mass frame.
   */
  const AngleDistribution& angle_distribution() const { return angle_; }

  /**
   * @breif Returns the Atomic Weight Ratio for the nuclide.
   */
  double awr() const { return awr_; }

  /**
   * @breif Returns the temperature for the nuclide.
   */
  double temperature() const { return awr_; }

  /**
   * @brief Returns the threshold for the application of the Target At Rest
   *        approximation.
   */
  double tar_threshold() const { return tar_threshold_; }

 private:
  struct Vector {
    double x, y, z;

    Vector(double x, double y, double z): x(x), y(y), z(z) {}

    Vector operator+(const Vector& v) const {
      return {x + v.x, y + v.y, z + v.z};
    }

    Vector operator-(const Vector& v) const {
      return {x - v.x, y - v.y, z - v.z};
    }

    Vector operator*(const double& c) const { return {x * c, y * c, z * c}; }

    Vector operator/(const double& c) const { return {x / c, y / c, z / c}; }

    double dot(const Vector& v) const { return x * v.x + y * v.y + z * v.z; }

    double magnitude() const;

    Vector rotate(double mu, double phi) const;
  };

  AngleDistribution angle_;
  double awr_;
  double temperature_;
  double tar_threshold_;

  Vector sample_target_velocity(double Ein, std::function<double()> rng) const;
};

}  // namespace pndl

#endif
