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
#ifndef PAPILLON_NDL_ANGLE_DISTRIBUTION_H
#define PAPILLON_NDL_ANGLE_DISTRIBUTION_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_law.hpp>
#include <functional>
#include <memory>

namespace pndl {

/**
 *  @brief Holds all of the angular distributions at all provided energies
 *         for a single reaction.
 */
class AngleDistribution {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param locb Index in the XSS array to start reading angular distribution.
   */
  AngleDistribution(const ACE& ace, int locb);

  /**
   * @param energy_grid Incoming energies which has angle laws.
   * @param laws Angle laws for each incoming energy.
   */
  AngleDistribution(const std::vector<double>& energy_grid,
                    const std::vector<std::shared_ptr<AngleLaw>>& laws);
  ~AngleDistribution() = default;

  /**
   * @brief Samples a scattering cosine for the given energy.
   * @param E_in Incident energy before scatter, in MeV.
   * @param rng Random number generator function.
   */
  double sample_angle(double E_in, std::function<double()> rng) const;

  /**
   * @brief Evaluates the PDF for having a scattering cosine of mu at incoming
   *        energy E_in.
   * @param E_in Incoming energy.
   * @param mu Scattering cosine.
   */
  double pdf(double E_in, double mu) const;

  /**
   * @brief Returns the number of energies/angular distributions stored.
   */
  std::size_t size() const { return energy_grid_.size(); }

  /**
   * @brief Reference to the vector of energy values (in MeV)
   *        which have an angular distribution.
   */
  const std::vector<double>& energy() const { return energy_grid_; }

  /**
   * @brief Gets the ith energy point (in MeV) which has an angular
   * distribution.
   * @param i Index in the energy grid.
   */
  double energy(std::size_t i) const { return energy_grid_[i]; }

  /**
   * @brief Gets a pointer to the angular distribution for the ith energy point.
   * @param i Index in the energy grid.
   */
  const AngleLaw& law(std::size_t i) const { return *laws_[i]; }

 private:
  std::vector<double> energy_grid_;
  std::vector<std::shared_ptr<AngleLaw>> laws_;
};

}  // namespace pndl

#endif
