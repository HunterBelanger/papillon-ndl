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
#include <PapillonNDL/delayed_family.hpp>
#include <PapillonNDL/equiprobable_energy_bins.hpp>
#include <PapillonNDL/evaporation.hpp>
#include <PapillonNDL/general_evaporation.hpp>
#include <PapillonNDL/maxwellian.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/tabular_energy.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <PapillonNDL/watt.hpp>
#include <iostream>

#include "constants.hpp"

namespace pndl {

DelayedFamily::DelayedFamily(const ACE& ace, std::size_t i, std::size_t g)
    : decay_constant_(), probability_(nullptr), energy_(nullptr) {
  // Read Decay Constant for family
  decay_constant_ = ace.xss(i);
  decay_constant_ *= SHAKE_TO_SEC;
  i++;

  // Read probability function
  uint32_t NR = ace.xss<uint32_t>(i);
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<double> energy = ace.xss(i + 2 + 2 * NR, NE);
  std::vector<double> y = ace.xss(i + 2 + 2 * NR + NE, NE);

  if (NR == 0 || NR == 1) {
    Interpolation interp = Interpolation::LinLin;
    if (NR == 1) interp = ace.xss<Interpolation>(i + 2);

    probability_ = std::make_shared<Tabulated1D>(interp, energy, y);
  } else {
    std::vector<uint32_t> breaks = ace.xss<uint32_t>(i + 1, NR);
    std::vector<Interpolation> interps = ace.xss<Interpolation>(i + 1 + NR, NR);

    probability_ = std::make_shared<Tabulated1D>(breaks, interps, energy, y);
  }

  // Get energy distribution location
  uint32_t locc = ace.xss<uint32_t>(static_cast<std::size_t>(ace.DNEDL()) + g - 1);
  std::size_t l = static_cast<std::size_t>(ace.DNED()) + locc - 1;

  // TODO currently ignore extra distribuitons, only read first one
  // Write a warning for now
  if (ace.xss<int>(l) != 0) {
    // there are other distributions
    std::cerr << "\n PapillonNDL WARNING : Delayed family " + std::to_string(g);
    std::cerr << " for ZAID " << ace.zaid() << " has multiple";
    std::cerr << " energy distributions.\n";
  }

  int law = ace.xss<int>(l + 1);
  uint32_t idat = ace.xss<uint32_t>(l + 2);
  std::size_t j = static_cast<std::size_t>(ace.DNED()) + idat - 1;

  if (law == 1) {  // Equiprobable Energy Bins
    energy_ = std::make_shared<EquiprobableEnergyBins>(ace, j);

  } else if (law == 4) {  // Tabular Energy
    energy_ = std::make_shared<TabularEnergy>(ace, j, ace.DNED());

  } else if (law == 5) {  // General Evaporation
    energy_ = std::make_shared<GeneralEvaporation>(ace, j);

  } else if (law == 7) {  // Maxwellian
    energy_ = std::make_shared<Maxwellian>(ace, j);

  } else if (law == 9) {  // Evaporation
    energy_ = std::make_shared<Evaporation>(ace, j);

  } else if (law == 11) {  // Watt
    energy_ = std::make_shared<Watt>(ace, j);

  } else {
    // Unknown or unsuported law
    std::string mssg = "Family " + std::to_string(g) +
                       " has unkown energy law " + std::to_string(law) + ".";
    throw PNDLException(mssg);
  }
}

}  // namespace pndl
