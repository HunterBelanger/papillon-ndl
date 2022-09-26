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
#ifndef PAPILLON_NDL_LEVEL_INELASTIC_SCATTER_H
#define PAPILLON_NDL_LEVEL_INELASTIC_SCATTER_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <memory>
#include <optional>

namespace pndl {

/**
 * @brief Energy distribution for inelastic scatter.
 */
class LevelInelasticScatter : public EnergyLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   */
  LevelInelasticScatter(const ACE& ace, std::size_t i);

  /**
   * @param Q Q-value of the reaction.
   * @param AWR Atomic Weight Ratio of nuclide.
   */
  LevelInelasticScatter(double Q, double AWR);
  ~LevelInelasticScatter() = default;

  double sample_energy(double E_in,
                       std::function<double()> rng) const override final;

  std::optional<double> pdf(double E_in, double E_out) const override final;

  /**
   * @brief Returns first parameter which is -(A+1)*Q/A.
   */
  double C1() const { return C1_; }

  /**
   * @brief Returns second parameter which is (A/(A+1))^2.
   */
  double C2() const { return C2_; }

 private:
  double C1_;
  double C2_;
};

}  // namespace pndl

#endif
