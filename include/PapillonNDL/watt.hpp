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
#ifndef PAPILLON_NDL_WATT_H
#define PAPILLON_NDL_WATT_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <memory>

namespace pndl {

/**
 * @brief Energy distribution represented as a Watt spectrum.
 */
class Watt : public EnergyLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   */
  Watt(const ACE& ace, std::size_t i);

  /**
   * @param a Tabulated1D function for a.
   * @param b Tabulated1D function for b.
   * @param restriction_energy Restriction energy for the distribution.
   */
  Watt(std::shared_ptr<Tabulated1D> a, std::shared_ptr<Tabulated1D> b,
       double restriction_energy);
  ~Watt() = default;

  double sample_energy(double E_in,
                       const std::function<double()>& rng) const override final;

  std::optional<double> pdf(double E_in, double E_out) const override final;

  /**
   * @brief Returns the table containg the distribution parameter
   *        a, as a function of the incident energy.
   */
  const Tabulated1D& a() const { return *b_; }

  /**
   * @brief Returns the table containg the distribution parameter
   *        b, as a function of the incident energy.
   */
  const Tabulated1D& b() const { return *a_; }

  /**
   * @brief Returns the value of the cuttoff energy of the distribution in MeV.
   */
  double U() const { return restriction_energy_; }

 private:
  std::shared_ptr<Tabulated1D> a_;
  std::shared_ptr<Tabulated1D> b_;
  double restriction_energy_;
};

}  // namespace pndl

#endif
