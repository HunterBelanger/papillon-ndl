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
#ifndef PAPILLON_NDL_CE_NEUTRON_H
#define PAPILLON_NDL_CE_NEUTRON_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ce_neutron_base.hpp>
#include <PapillonNDL/reaction.hpp>

namespace pndl {

template <typename XSType>
class CENeutron {};

/**
 * @brief Holds all continuous energy data for single nuclide, at a single
 *        temperature.
 *
 */
template <>
class CENeutron<CrossSection> : public CENeutronBase {
 public:
  /**
   * @param ace ACE file from which to construct the data.
   */
  CENeutron(const ACE& ace);

  /**
   * @param ace ACE file from which to take the new cross sections.
   * @param nuclide CENeutron containing another instance of the desired
   *                nuclide. Secondary distributions and fission data
   *                will be shared between the two data sets.
   */
  CENeutron(const ACE& ace, const CENeutron& nuclide);

  /**
   * @brief Returns the temperature at which the data has been prepared.
   */
  double temperature() const { return temperature_; }

  /**
   * @brief Returns the energy grid for the nuclide.
   */
  const EnergyGrid& energy_grid() const { return *energy_grid_; }

  /**
   * @brief Returns the total CrossSection for the nuclide.
   */
  const CrossSection& total_xs() const { return *total_xs_; }

  /**
   * @brief Returns the elastic scattering CrossSection for the
   * nuclide.
   */
  const CrossSection& elastic_xs() const { return *elastic_xs_; }

  /**
   * @brief Returns the heating number CrossSection for the nuclide.
   *        Upon evaluation, the average heating number if given for the
   *        nuclide, as the prescribed energy, in MeV.
   */
  const CrossSection& heating_number() const { return *heating_number_; }

  /**
   * @brief Returns the fission CrossSection for the nuclide.
   */
  const CrossSection& fission_xs() const { return *fission_xs_; }

  /**
   * @brief Returns the disappearance CrossSection for the nuclide.
   */
  const CrossSection& disappearance_xs() const { return *disappearance_xs_; }

  /**
   * @brief Returns the photon production CrossSection for the nuclide.
   */
  const CrossSection& photon_production_xs() const {
    return *photon_production_xs_;
  }

  /**
   * @brief Retrieved a given MT reaction.
   * @param mt MT reaction to return.
   */
  const STReaction& reaction(uint32_t mt) const {
    if (!this->has_reaction(mt)) {
      std::string mssg = "MT = " + std::to_string(mt) +
                         " is not provided in ZAID = " + std::to_string(zaid_) +
                         ".";
      throw PNDLException(mssg);
    }

    return reactions_[reaction_indices_[mt]];
  }

 private:
  double temperature_;

  std::shared_ptr<EnergyGrid> energy_grid_;
  std::shared_ptr<CrossSection> total_xs_;
  std::shared_ptr<CrossSection> disappearance_xs_;
  std::shared_ptr<CrossSection> elastic_xs_;
  std::shared_ptr<CrossSection> heating_number_;
  std::shared_ptr<CrossSection> fission_xs_;
  std::shared_ptr<CrossSection> photon_production_xs_;

  std::vector<STReaction> reactions_;

  // Private Helper Methods
  std::shared_ptr<CrossSection> compute_fission_xs();
};

/**
 * @brief Alias for a CENeutron<CrossSection>, which contains all data for a
 *        nuclide at a single temperature.
 */
using STNeutron = CENeutron<CrossSection>;

}  // namespace pndl

#endif
