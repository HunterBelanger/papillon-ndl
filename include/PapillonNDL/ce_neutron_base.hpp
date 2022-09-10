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
#ifndef PAPILLON_NDL_CE_NEUTRON_BASE_H
#define PAPILLON_NDL_CE_NEUTRON_BASE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_distribution.hpp>
#include <PapillonNDL/delayed_group.hpp>
#include <PapillonNDL/nuclide.hpp>
#include <PapillonNDL/zaid.hpp>
#include <array>
#include <memory>

namespace pndl {

/**
 * @brief Holds all non temperature-dependent, continuous energy data for a
 *        single nuclide. This is mainly ZAID, AWR, delayed fission info,
 *        and nu.
 */
class CENeutronBase {
 public:
  CENeutronBase(const CENeutronBase&) = default;

  /**
   * @brief Returns the nuclide ZAID.
   */
  const ZAID& zaid() const { return zaid_; }

  /**
   * @brief Returns the nuclide Atomic Weight Ratio.
   */
  double awr() const { return awr_; }

  /**
   * @brief Returns true if the nuclide is fissile, and false otherwise.
   */
  bool fissile() const { return fissile_; }

  /**
   * @brief Returns a list of all scattering and absorption MT reactions present
   *        for the nuclide (other than elastic).
   */
  const std::vector<uint32_t>& mt_list() const { return mt_list_; }

  /**
   * @brief Checks to see if a nucldie has a given scattering or absorption
   *        reaction.
   * @param mt MT reaction to search for.
   */
  bool has_reaction(uint32_t mt) const {
    return (mt > 891 || reaction_indices_[mt] < 0) ? false : true;
  }

 protected:
  ZAID zaid_;
  double awr_;
  bool fissile_;

  std::vector<uint32_t> mt_list_;
  std::array<int32_t, 892> reaction_indices_;

  /**
   * @param ace ACE file from which to construct the data.
   */
  CENeutronBase(const ACE& ace);
};

}  // namespace pndl

#endif
