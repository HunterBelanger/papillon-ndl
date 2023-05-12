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
#include <PapillonNDL/continuous_energy_discrete_cosines.hpp>
#include <PapillonNDL/discrete_cosines_energies.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/st_incoherent_inelastic.hpp>

namespace pndl {

STIncoherentInelastic::STIncoherentInelastic(const ACE& ace,
                                             bool unit_based_interpolation)
    : xs_(nullptr), angle_energy_(nullptr) {
  // Read the XS
  try {
    std::size_t S = static_cast<std::size_t>(ace.jxs(0) - 1);
    uint32_t Ne = ace.xss<uint32_t>(S);  // Number of grid points
    std::vector<double> energy = ace.xss(S + 1, Ne);
    std::vector<double> xs = ace.xss(S + 1 + Ne, Ne);
    xs_ = std::make_shared<Tabulated1D>(Interpolation::LinLin, energy, xs);
  } catch (PNDLException& err) {
    std::string mssg = "Could not construct cross section.";
    err.add_to_exception(mssg);
    throw err;
  }

  // Read the angle-energy distribution
  try {
    int32_t nxs_7 = ace.nxs(6);
    if (nxs_7 == 0 || nxs_7 == 1) {
      angle_energy_ = std::make_shared<DiscreteCosinesEnergies>(ace);
    } else if (nxs_7 == 2) {
      angle_energy_ = std::make_shared<ContinuousEnergyDiscreteCosines>(
          ace, unit_based_interpolation);
    } else {
      std::string mssg =
          "Unknown distribution type. Make sure this is a valid thermal "
          "scattering law ACE file.";
      throw PNDLException(mssg);
    }
  } catch (PNDLException& err) {
    std::string mssg = "Could not construct AngleEnergy distribution.";
    err.add_to_exception(mssg);
    throw err;
  }
}
}  // namespace pndl
