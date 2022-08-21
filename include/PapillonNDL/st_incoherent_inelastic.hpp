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
#ifndef PAPILLON_NDL_ST_INCOHERENT_INELASTIC_H
#define PAPILLON_NDL_ST_INCOHERENT_INELASTIC_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/tabulated_1d.hpp>

namespace pndl {

/**
 * @brief Holds the Incoherent Inelastic scattering data for a single nuclide
 *        at a single temperature.
 */
class STIncoherentInelastic {
 public:
  /**
   * @param ace ACE file which contains thermal scattering law.
   * @param unit_based_interpolation If false (default value), the distribution
   *        will be sampled without using unit-based interpolation, which is
   *        the method used by MCNP, Serpent, and OpenMC. If set to true, unit
   *        based interpolation will be applied to the sampling of the energy.
   */
  STIncoherentInelastic(const ACE& ace, bool unit_based_interpolation = false);
  ~STIncoherentInelastic() = default;

  /**
   * @brief Returns the cross section function.
   */
  const Tabulated1D& xs() const { return *xs_; }

  /**
   * @brief Evaluates the incoherent inelastic scattering cross section
   *        at energy E.
   * @param E Incident energy at which to evaluate the cross section in MeV.
   */
  double xs(double E) const { return xs_->evaluate(E); }

  /**
   * @brief Retruns the maximum energy value which is tabulated for the
   *        cross section.
   */
  double max_energy() const { return this->xs_->max_x(); }

  /**
   * @brief Sample the angle-energy distribution.
   * @param E_in Incident energy in MeV.
   * @param rng Random number generation function.
   */
  AngleEnergyPacket sample_angle_energy(double E_in,
                                        std::function<double()> rng) const {
    return angle_energy_->sample_angle_energy(E_in, rng);
  }

  /**
   * @brief Returns the AngleEnergy distribution.
   */
  const AngleEnergy& distribution() const { return *angle_energy_; }

 private:
  std::shared_ptr<Tabulated1D> xs_;
  std::shared_ptr<AngleEnergy> angle_energy_;
};

}  // namespace pndl

#endif
