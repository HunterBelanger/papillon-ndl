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
#ifndef PAPILLON_NDL_TABULAR_ENERGY_H
#define PAPILLON_NDL_TABULAR_ENERGY_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <PapillonNDL/pctable.hpp>

namespace pndl {

/**
 * @brief Energy distribution represented as a tabulated PDF and CDF.
 */
class TabularEnergy : public EnergyLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   * @param JED Index offset to find the PDF and CDF tables.
   */
  TabularEnergy(const ACE& ace, std::size_t i, std::size_t JED);

  /**
   * @param incoming_energy Incoming energy grid.
   * @param tables PCTables for each point in the incoming energy grid.
   */
  TabularEnergy(const std::vector<double>& incoming_energy,
                const std::vector<PCTable>& tables);

  ~TabularEnergy() = default;

  double sample_energy(double E_in,
                       std::function<double()> rng) const override final;

  std::optional<double> pdf(double E_in, double E_out) const override final;

  /**
   * @brief Reterns the incoming energy points in MeV for which a PCTable
   *        is stored.
   */
  const std::vector<double>& incoming_energy() const {
    return incoming_energy_;
  }

  /**
   * @brief Returns the ith PDF/CDF, corresponding to the ith incoming energy.
   * @param i Index in the incoming energy grid.
   */
  const PCTable& table(std::size_t i) const { return tables_[i]; }

  /**
   * @brief Returns the number of incoming energy points.
   */
  std::size_t size() const { return incoming_energy_.size(); }

 private:
  std::vector<double> incoming_energy_;
  std::vector<PCTable> tables_;
};

}  // namespace pndl

#endif
