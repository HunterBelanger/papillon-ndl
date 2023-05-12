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
#ifndef PAPILLON_NDL_ST_TSL_REACTION_H
#define PAPILLON_NDL_ST_TSL_REACTION_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/angle_energy.hpp>

namespace pndl {

/**
 * @brief Interface to represent any thermal neutron scattering reaction at a
 *        single temperature.
 */
class STTSLReaction : public AngleEnergy {
 public:
  /**
   * @brief Evaluates the cross section for the thermal neutron scattering
   *        reaction.
   * @param E Incident energy at which to evaluate the cross section in MeV.
   */
  virtual double xs(double E) const = 0;
};

}  // namespace pndl

#endif
