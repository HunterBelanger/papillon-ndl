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
#ifndef PAPILLON_NDL_ST_NEUTRON_H
#define PAPILLON_NDL_ST_NEUTRON_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/elastic.hpp>
#include <PapillonNDL/fission.hpp>
#include <PapillonNDL/reaction.hpp>
#include <PapillonNDL/urr_ptables.hpp>
#include <PapillonNDL/xs_packet.hpp>
#include <memory>

namespace pndl {

/**
 * @brief Holds all continuous energy data for single nuclide, at a single
 *        temperature.
 *
 */
class STNeutron {
 public:
  /**
   * @param ace ACE file from which to construct the data.
   */
  STNeutron(const ACE& ace);

  /**
   * @param ace ACE file from which to take the new cross sections.
   * @param nuclide CENeutron containing another instance of the desired
   *                nuclide. Secondary distributions and fission data
   *                will be shared between the two data sets.
   */
  STNeutron(const ACE& ace, const STNeutron& nuclide);

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
    if (mt == 18 || mt == 19 || mt == 20 || mt == 21 || mt == 38) {
      return fission_->has_reaction(mt);
    } else {
      return (mt > 891 || reaction_indices_[mt] < 0) ? false : true;
    }
  }

  /**
   * @brief Retrieved a given MT reaction.
   * @param mt MT reaction to return.
   */
  const STReaction& reaction(uint32_t mt) const {
    if (!this->has_reaction(mt)) {
      std::string mssg =
          "MT = " + std::to_string(mt) +
          " is not provided in ZAID = " + std::to_string(zaid_.zaid()) + ".";
      throw PNDLException(mssg);
    }

    if (mt == 18 || mt == 19 || mt == 20 || mt == 21 || mt == 38) {
      return fission_->reaction(mt);
    } else {
      return reactions_[static_cast<std::size_t>(reaction_indices_[mt])];
    }
  }

  /**
   * @brief Returns a reference to the URRPTables instance.
   */
  const URRPTables& urr_ptables() const { return *urr_ptables_; }

  /**
   * @brief Returns a reference to the Elastic instance which contains the
   *        AngleEnergy distribution for elastic scattering.
   */
  const Elastic& elastic() const { return *elastic_; }

  /**
   * @brief Returns a modifiable reference to the Elastic instance which
   *        contains the AngleEnergy distribution for elastic scattering.
   */
  Elastic& elastic() { return *elastic_; }

  /**
   * @brief Returns and reference to the Fission instance which contains all
   *        fission information.
   */
  const Fission& fission() { return *fission_; }

  /**
   * @brief Evaluates the important nuclide cross sections at a given energy,
   *        with the grid point already provided.
   * @param E Energy to evaluate the cross section at.
   * @param i Index of the points for interpolation in the frame of the energy
   *          grid.
   */
  XSPacket evaluate_xs(double Ein, std::size_t i) const {
    XSPacket xs;
    xs.total = total_xs_->evaluate(Ein, i);
    xs.elastic = elastic_xs_->evaluate(Ein, i);
    xs.fission = fission_xs_->evaluate(Ein, i);
    xs.absorption = disappearance_xs_->evaluate(Ein, i) + xs.fission;
    xs.heating = heating_number_->evaluate(Ein, i);
    xs.inelastic = xs.total - xs.elastic - xs.absorption;

    if (xs.inelastic < 0.) xs.inelastic = 0.;

    if (this->has_reaction(102)) {
      xs.capture = this->reaction(102).xs()(Ein, i);
    } else {
      xs.capture = 0.;
    }

    return xs;
  }

  /**
   * @brief Evaluates the important nuclide cross sections at a given energy.
   * @param E Energy to evaluate the cross section at.
   */
  XSPacket evaluate_xs(double Ein) const {
    std::size_t i = energy_grid_->get_lower_index(Ein);
    return this->evaluate_xs(Ein, i);
  }

 private:
  ZAID zaid_;
  double awr_;
  bool fissile_;
  double temperature_;

  std::shared_ptr<EnergyGrid> energy_grid_;
  std::shared_ptr<CrossSection> total_xs_;
  std::shared_ptr<CrossSection> disappearance_xs_;
  std::shared_ptr<CrossSection> elastic_xs_;
  std::shared_ptr<CrossSection> heating_number_;
  std::shared_ptr<CrossSection> fission_xs_;
  std::shared_ptr<CrossSection> photon_production_xs_;

  std::shared_ptr<Elastic> elastic_;
  std::shared_ptr<Fission> fission_;
  std::shared_ptr<URRPTables> urr_ptables_;

  std::vector<uint32_t> mt_list_;
  std::array<int32_t, 892> reaction_indices_;
  std::vector<STReaction> reactions_;

  // Private Helper Methods
  std::shared_ptr<CrossSection> compute_fission_xs();
};

}  // namespace pndl

#endif
