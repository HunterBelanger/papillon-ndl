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
#ifndef PAPILLON_NDL_CONTINUOUS_ENERGY_DISCRETE_COSINES_H
#define PAPILLON_NDL_CONTINUOUS_ENERGY_DISCRETE_COSINES_H

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>

namespace pndl {

/**
 * @brief Class which represents continuous energy distributions with
 *        discrete cosines for incoherent inelastic scattering.
 */
class ContinuousEnergyDiscreteCosines : public AngleEnergy {
 public:
  /**
   * @param ace ACE file which contains thermal scattering law.
   * @param unit_based_interpolation If false (default value), the distribution
   *        will be sampled without using unit-based interpolation, which is
   *        the method used by MCNP, Serpent, and OpenMC. If set to true, unit
   *        based interpolation will be applied to the sampling of the energy.
   */
  ContinuousEnergyDiscreteCosines(const ACE& ace,
                                  bool unit_based_interpolation = false);

  /**
   * @brief Struct which contains the outgoing energy distribution and
   *        the discrete scattering cosiens for a single incident energy.
   */
  struct CEDCTable {
    std::vector<double> energy; /**< Outgoing energy points */
    std::vector<double> pdf;    /**< PDF for the outgoing energy */
    std::vector<double> cdf;    /**< CDF for the outgoing energy */
    std::vector<std::vector<double>>
        cosines; /**< Discrete scattering cosines for each outgoing energy */

    /**
     * @brief Samples and outgoing energy from the distributions while
     *        also setting the value j to be the lower bound index
     *        for sampling the scattering cosine latter on.
     * @param xi Random variable on the unit interval [0,1).
     * @param j Index which will be set to latter locate the proper
     *          distribution for the scattering cosine.
     */
    double sample_energy(double xi, std::size_t& j) const;
  };

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  std::optional<double> angle_pdf(double E_in, double mu) const override final;

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override final;

  /**
   * @brief Returns vector to the incoming energy grid.
   */
  const std::vector<double>& incoming_energy() const {
    return incoming_energy_;
  }

  /**
   * @brief Returns the number of incoming energy points.
   */
  std::size_t size() const { return incoming_energy_.size(); }

  /**
   * @brief Returns the vector of all CEDCTables.
   */
  const std::vector<CEDCTable> tables() const { return tables_; }

  /**
   * @brief Returns a CEDCTable which contains the distributions
   *        for the ith incoming energy.
   * @param i Index to the incoming energy.
   */
  const CEDCTable& table(std::size_t i) const { return tables_[i]; }

  /**
   * @brief Returns true if the distribution used unit-based interpolation
   *        in sampling the scattering energy and angle, and false otherwise.
   */
  bool unit_based_interpolation() const { return unit_based_interpolation_; }

 private:
  std::vector<double> incoming_energy_;
  std::vector<CEDCTable> tables_;
  uint32_t Nmu;
  bool unit_based_interpolation_;

  AngleEnergyPacket sample_with_unit_based_interpolation(
      double E_in, std::function<double()> rng) const;
  AngleEnergyPacket sample_without_unit_based_interpolation(
      double E_in, std::function<double()> rng) const;
};

}  // namespace pndl

#endif
