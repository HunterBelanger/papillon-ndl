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
#ifndef PAPILLON_NDL_FISSION_H
#define PAPILLON_NDL_FISSION_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/delayed_family.hpp>
#include <PapillonNDL/energy_grid.hpp>
#include <PapillonNDL/function_1d.hpp>
#include <PapillonNDL/reaction.hpp>
#include <PapillonNDL/zaid.hpp>
#include <memory>
#include <vector>

namespace pndl {

/**
 * @brief This class contains all of the information for the nuclide which is
 * related to
 */
class Fission {
 public:
  /**
   * @param ace ACE file from which to construct the data.
   * @param energy_grid Pointer to the EnergyGrid for the nuclide.
   */
  Fission(const ACE& ace, std::shared_ptr<EnergyGrid> energy_grid);

  /**
   * @param ace ACE file from which to construct the data.
   * @param energy_grid Pointer to the EnergyGrid for the nuclide.
   * @param fission Fission object to take distributions from.
   */
  Fission(const ACE& ace, std::shared_ptr<EnergyGrid> energy_grid,
          const Fission& fission);

  /**
   * @brief Returns the function for total nu.
   */
  const Function1D& nu_total() const { return *nu_total_; }

  /**
   * @brief Returns the function for prompt nu.
   */
  const Function1D& nu_prompt() const { return *nu_prompt_; }

  /**
   * @brief Returns the function for delayed nu.
   */
  const Function1D& nu_delayed() const { return *nu_delayed_; }

  /**
   * @brief Returns the prompt spectrum for fission neutrons.
   */
  const AngleEnergy& prompt_spectrum() const { return *prompt_spectrum_; }

  /**
   * @brief Returns the number of delayed neutron families.
   */
  std::size_t n_delayed_families() const { return delayed_families_.size(); }

  /**
   * @brief Returns the ith delayed family data.
   * @param i Index of the delayed family.
   */
  const DelayedFamily& delayed_family(std::size_t i) const {
    return delayed_families_[i];
  }

  /**
   * @brief Returns a list of fission reactions present.
   */
  const std::vector<uint32_t>& mt_list() const { return mt_list_; }

  /**
   * @brief Checks to see if a given fission reaction is present. The only MT
   *        values which could possibly be present are 18, 19, 20, 21, and 38.
   *        All other MTs are guaranteed to return false.
   * @param mt MT fission reaction to search for.
   */
  bool has_reaction(uint32_t mt) const {
    if (mt == 18 && mt18_)
      return true;
    else if (mt == 19 && mt19_)
      return true;
    else if (mt == 20 && mt20_)
      return true;
    else if (mt == 21 && mt21_)
      return true;
    else if (mt == 38 && mt38_)
      return true;
    else
      return false;
  }

  /**
   * @brief Retrieves a given MT fission reaction. The only MT values which
   *        could possibly be present are 18, 19, 20, 21, and 38 (check with
   *        has_reaction). All other MTs are guaranteed to throw a
   *        PNDLException.
   * @param mt MT fission reaction to return.
   */
  const STReaction& reaction(uint32_t mt) const {
    if (!this->has_reaction(mt)) {
      std::string mssg =
          "MT = " + std::to_string(mt) +
          " is not provided in ZAID = " + std::to_string(zaid_.zaid()) + ".";
      throw PNDLException(mssg);
    }

    if (mt == 18)
      return *mt18_;
    else if (mt == 19)
      return *mt19_;
    else if (mt == 20)
      return *mt20_;
    else if (mt == 21)
      return *mt21_;
    else
      return *mt38_;
  }

 private:
  ZAID zaid_;
  std::shared_ptr<Function1D> nu_total_;
  std::shared_ptr<Function1D> nu_prompt_;
  std::shared_ptr<Function1D> nu_delayed_;
  std::shared_ptr<STReaction> mt18_;
  std::shared_ptr<STReaction> mt19_;
  std::shared_ptr<STReaction> mt20_;
  std::shared_ptr<STReaction> mt21_;
  std::shared_ptr<STReaction> mt38_;
  std::shared_ptr<const AngleEnergy> prompt_spectrum_;
  std::vector<DelayedFamily> delayed_families_;
  std::vector<uint32_t> mt_list_;

  // Private helper methods
  std::shared_ptr<Function1D> read_nu(const ACE& ace, std::size_t i);
  std::shared_ptr<Function1D> read_polynomial_nu(const ACE& ace, std::size_t i);
  std::shared_ptr<Function1D> read_tabular_nu(const ACE& ace, std::size_t i);
};

}  // namespace pndl

#endif
