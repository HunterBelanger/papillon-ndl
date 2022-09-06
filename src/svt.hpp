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
#ifndef PAPILLON_NDL_SVT_H
#define PAPILLON_NDL_SVT_H

#include <PapillonNDL/cross_section.hpp>
#include <functional>

#include "vector.hpp"

namespace pndl {

/**
 * @brief Samples the velocity of a target nuclide from a Maxwelliam spectrum,
 *        while assuming that the elastic scatting cross section is constant in
 *        the vicinity of Ein. It is always assumed that the direction of the
 *        incident neutron is (0,0,1).
 * @param Ein Incident energy of the neutron in MeV.
 * @param kT Temperature of the "free-gas" in MeV.
 * @param awr Atomic weight ratio of the nuclide.
 * @param rng Random number generator function.
 */
Vector sample_target_velocity(const double& Ein, const double& kT,
                              const double& awr,
                              const std::function<double()>& rng);

/**
 * @brief Samples the velocity of a target nuclide from a Maxwelliam spectrum,
 *        while assuming that the elastic scatting cross section is constant in
 *        the vicinity of Ein. It is always assumed that the direction of the
 *        incident neutron is (0,0,1).
 * @param Ein Incident energy of the neutron in MeV.
 * @param xs Elastic scattering CrossSection at 0 Kelvin.
 * @param kT Temperature of the "free-gas" in MeV.
 * @param awr Atomic weight ratio of the nuclide.
 * @param rng Random number generator function.
 */
Vector sample_target_velocity_dbrc(const double& Ein, const CrossSection& xs,
                                   const double& kT, const double& awr,
                                   const std::function<double()>& rng);

}  // namespace pndl

#endif
