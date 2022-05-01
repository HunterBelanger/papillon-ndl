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
#ifndef PAPILLON_NDL_REACTION_H
#define PAPILLON_NDL_REACTION_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/reaction_base.hpp>
#include <memory>

namespace pndl {

template <typename XSType>
class Reaction {};

/**
 * @brief Holds all information for a reaction at a single temperature.
 */
template <>
class Reaction<CrossSection> : public ReactionBase {
 public:
  /**
   * @param ace ACE file to take reaction from.
   * @param indx Reaction index in the MT array.
   * @param egrid EnergyGrid for the nuclide.
   */
  Reaction(const ACE& ace, std::size_t indx, const EnergyGrid& egrid);

  /**
   * @param ace ACE file to take cross section from.
   * @param indx Reaction index in the MT array.
   * @param egrid EnergyGrid for the nuclide.
   * @param reac Reaction object to take distributions from.
   */
  Reaction(const ACE& ace, std::size_t indx, const EnergyGrid& egrid,
           const Reaction& reac);

  /**
   * @brief Returns the CrossSection for the reaction.
   */
  const CrossSection& xs() const { return *xs_; }

 private:
  std::shared_ptr<CrossSection> xs_;
};

/**
 * @brief Alias for Reaction<CrossSection>, which contains all data for a
 *        reaction MT at a single temperature.
 */
using STReaction = Reaction<CrossSection>;

}  // namespace pndl

#endif
