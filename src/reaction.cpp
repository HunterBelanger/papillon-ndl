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
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/reaction.hpp>
#include <memory>

/**
 * @file
 * @author Hunter Belanger
 */

namespace pndl {

Reaction<CrossSection>::Reaction(const ACE& ace, std::size_t indx,
                                 std::shared_ptr<EnergyGrid> egrid)
    : ReactionBase(ace, indx), xs_(nullptr) {
  try {
    uint32_t loca = ace.xss<uint32_t>(ace.LSIG() + indx);
    xs_ = std::make_shared<CrossSection>(ace, ace.SIG() + loca - 1, egrid);
    threshold_ = xs_->energy(0);
  } catch (PNDLException& error) {
    std::string mssg = "Could not create cross section for MT = " +
                       std::to_string(this->mt()) + ".";
    error.add_to_exception(mssg);
    throw error;
  }
}

Reaction<CrossSection>::Reaction(const ACE& ace, std::size_t indx,
                                 std::shared_ptr<EnergyGrid> egrid,
                                 const Reaction& reac)
    : ReactionBase(reac), xs_(nullptr) {
  // make sure the MT values agree
  if (this->mt() != reac.mt()) {
    std::string mssg = "MT from ACE file doesn't match MT from reaction.";
    throw PNDLException(mssg);
  }

  // Get XS from new ACE
  try {
    uint32_t loca = ace.xss<uint32_t>(ace.LSIG() + indx);
    xs_ = std::make_shared<CrossSection>(ace, ace.SIG() + loca - 1, egrid);
    threshold_ = xs_->energy(0);
  } catch (PNDLException& error) {
    std::string mssg =
        "Could not create cross section for MT = " + std::to_string(mt_) + ".";
    error.add_to_exception(mssg);
    throw error;
  }
}

Reaction<CrossSection>::Reaction(
    const CrossSection& xs, uint32_t mt, double q, double awr, double threshold,
    std::shared_ptr<Function1D> yield,
    std::shared_ptr<AngleEnergy> neutron_distribution) try
    : ReactionBase(mt, q, awr, threshold, yield, neutron_distribution),
      xs_(nullptr) {
  xs_ = std::make_shared<CrossSection>(xs);
} catch (PNDLException& err) {
  std::string mssg = "Could not create ReactionBase.";
  err.add_to_exception(mssg);
  throw err;
}

}  // namespace pndl
