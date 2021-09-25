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
#ifndef PAPILLON_NDL_ABSORPTION_H
#define PAPILLON_NDL_ABSORPTION_H

#include <PapillonNDL/angle_energy.hpp>
#include <cstddef>
#include <optional>

/**
 * @file
 * @author Hunter Belanger
 */

namespace pndl {

/**
 * @brief A distribution to represent absorption.
 */
class Absorption : public AngleEnergy {
 public:
  Absorption() {}

  AngleEnergyPacket sample_angle_energy(
      double /*E_in*/, std::function<double()> /*rng*/) const override final {
    return {1., 0.};
  }

  std::optional<double> angle_pdf(double /*E_in*/,
                                  double /*mu*/) const override final {
    return std::nullopt;
  }

  std::optional<double> pdf(double /*E_in*/, double /*mu*/,
                            double /*E_out*/) const override final {
    return std::nullopt;
  }
};

}  // namespace pndl

#endif
