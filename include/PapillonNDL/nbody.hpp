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
#ifndef PAPILLON_NDL_NBODY_H
#define PAPILLON_NDL_NBODY_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>

namespace pndl {

/**
 * @brief Implements product Angle-Energy law which follows an N-body phase
 *        space distribution.
 */
class NBody : public AngleEnergy {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   * @param iQ Q-value for the reaction.
   */
  NBody(const ACE& ace, std::size_t i, double iQ);

  /**
   * @param n Number of particles (3, 4, or 5).
   * @param Ap Total mass ratio for the n particles.
   * @param AWR Atomic Weight Ratio of the nuclide.
   * @param Q Q-value for the reaction.
   */
  NBody(uint16_t n, double Ap, double AWR, double Q);
  ~NBody() = default;

  AngleEnergyPacket sample_angle_energy(
      double E_in, const std::function<double()>& rng) const override final;

  std::optional<double> angle_pdf(double E_in, double mu) const override final;

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override final;

  /**
   * @brief Returns the number of bodies.
   */
  uint32_t n() const { return n_; }

  /**
   * @brief Returns the total AWR for all of the particles.
   */
  double Ap() const { return Ap_; }

  /**
   * @brief Returns the AWR of the nuclide in question.
   */
  double A() const { return A_; }

  /**
   * @brief Returns the Q-value for the reaction in MeV.
   */
  double Q() const { return Q_; }

 private:
  uint32_t n_;
  double Ap_;
  double A_;
  double Q_;

  double maxwellian_spectrum(const std::function<double()>& rng) const;
};

}  // namespace pndl

#endif
