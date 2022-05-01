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
#include <PapillonNDL/s1_cross_section.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <memory>
#include <sstream>
#include <string>

namespace pndl {

S1CrossSection::S1CrossSection(const ACE& ace, std::size_t i,
                           const EnergyGrid& E_grid, bool get_index)
    try: xs_(ace, i, E_grid, get_index), min_temp_(ace.temperature()), awr_(ace.awr()) {
  if (min_temp_ < 0.) {
    std::stringstream mssg;
    mssg << "Temperature must be greater than or equal to zero.";
    mssg << " ACE file provided temperature = " << awr_ << ".";
    throw PNDLException(mssg.str());
  }

  if (awr_ <= 0.) {
    std::stringstream mssg;
    mssg << "AWR from provided ACE file is less than or equal to zero.";
    mssg << " Provided AWR = " << awr_ << ".";
    throw PNDLException(mssg.str());
  }
} catch (PNDLException& error) {
  error.add_to_exception("Could not construct S1CrossSection.");
  throw error;
}

S1CrossSection::S1CrossSection(const std::vector<double>& xs,
                           const EnergyGrid& E_grid,
                           std::size_t index,
                           double temperature,
                           double awr)
    try : xs_(xs,E_grid,index),
      min_temp_(temperature),
      awr_(awr) {
  if (awr_ <= 0.) {
    std::stringstream mssg;
    mssg << "AWR must be greater than zero.";
    mssg << " ACE file provided AWR = " << awr_ << ".";
    throw PNDLException(mssg.str());
  }

  if (min_temp_ < 0.) {
    std::stringstream mssg;
    mssg << "Temperature must be greater than or equal to zero.";
    mssg << " ACE file provided temperature = " << awr_ << ".";
    throw PNDLException(mssg.str());
  }
} catch (PNDLException& error) {
  error.add_to_exception("Could not construct S1CrossSection.");
  throw error;
}

S1CrossSection::S1CrossSection(double xs, const EnergyGrid& E_grid)
    try: xs_(xs, E_grid), min_temp_(0.), awr_(0.) {
    } catch (PNDLException& error) {
      error.add_to_exception("Could not construct S1CrossSection.");
      throw error;
    }

S1CrossSection::S1CrossSection(const CrossSection& xs, double temperature, double awr)
  try: xs_(xs), min_temp_(temperature), awr_(awr) {
  if (awr_ <= 0.) {
    std::stringstream mssg;
    mssg << "AWR must be greater than zero.";
    mssg << " ACE file provided AWR = " << awr_ << ".";
    throw PNDLException(mssg.str());
  }

  if (min_temp_ < 0.) {
    std::stringstream mssg;
    mssg << "Temperature must be greater than or equal to zero.";
    mssg << " ACE file provided temperature = " << awr_ << ".";
    throw PNDLException(mssg.str());
  }
} catch (PNDLException& error) {
  error.add_to_exception("Could not construct S1CrossSection.");
  throw error;
}

}  // namespace pndl
