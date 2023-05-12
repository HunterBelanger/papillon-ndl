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
#ifndef PAPILLON_NDL_ISOTROPIC_H
#define PAPILLON_NDL_ISOTROPIC_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/angle_law.hpp>
#include <cmath>
#include <functional>

namespace pndl {

/**
 * @brief Isotropic angular distribution.
 */
class Isotropic : public AngleLaw {
 public:
  Isotropic() {}
  ~Isotropic() = default;

  double sample_mu(const std::function<double()>& rng) const override final;

  double pdf(double mu) const override final;
};

}  // namespace pndl

#endif
