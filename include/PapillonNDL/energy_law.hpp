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
#ifndef PAPILLON_NDL_ENERGY_LAW_H
#define PAPILLON_NDL_ENERGY_LAW_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <functional>
#include <memory>
#include <optional>

namespace pndl {

/**
 * @brief Interface to represent uncorrelated energy distributions.
 */
class EnergyLaw : public std::enable_shared_from_this<EnergyLaw> {
 public:
  virtual ~EnergyLaw() = default;

  /**
   * @brief Samples an energy (in MeV) from the distribution.
   * @param E_in Incident energy in MeV.
   * @param rng Random number generation function.
   */
  virtual double sample_energy(double E_in,
                               const std::function<double()>& rng) const = 0;

  /**
   * @brief Samples the PDF for the energy transfer from E_in to E_out where
   *        E_in is provided in the lab frame, and E_out is provided in the
   *        frame of the reaction data. Returned as an std::optional<double>,
   *        as some EnergyLaw types have no defined PDF.
   * @param E_in Incoming energy.
   * @param E_out Outgoing energy.
   */
  virtual std::optional<double> pdf(double E_in, double E_out) const = 0;
};

}  // namespace pndl

#endif
