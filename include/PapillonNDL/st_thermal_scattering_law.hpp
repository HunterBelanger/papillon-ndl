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
#ifndef PAPILLON_NDL_ST_THERMAL_SCATTERING_LAW_H
#define PAPILLON_NDL_ST_THERMAL_SCATTERING_LAW_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/st_incoherent_inelastic.hpp>
#include <PapillonNDL/st_tsl_reaction.hpp>
#include <PapillonNDL/zaid.hpp>

namespace pndl {

/**
 * @brief Class to hold all thermal scattering data for a single nuclide
 *        at at single temperature.
 */
class STThermalScatteringLaw {
 public:
  /**
   * @param ace ACE file which contains thermal scattering law.
   * @param unit_based_interpolation If false (default value) and the incoherent
   *        inelastic scattering distribution is continuous in energy, unit
   *        based interpolation will not be applied. This is method used by
   *        MCNP, Serpent, and OpenMC, so we have made it our default. If set to
   *        true, unit based interpolation will be used.
   */
  STThermalScatteringLaw(const ACE& ace, bool unit_based_interpolation = false);
  ~STThermalScatteringLaw() = default;

  /**
   * @brief Returns the nuclide ZAID.
   */
  const ZAID& zaid() const { return zaid_; }

  /**
   * @brief Returns the nuclide Atomic Weight Ratio.
   */
  double awr() const { return awr_; }

  /**
   * @brief Returns the temperature at which the data has been prepared.
   */
  double temperature() const { return temperature_; }

  /**
   * @brief Returns the maximum energy for the incoherent inelastic scattering
   *        reaction. This value is typtically used as the cutoff for using
   *        Sab tables in Monte Carlo codes.
   */
  double max_energy() const { return incoherent_inelastic_->max_energy(); }

  /**
   * @brief Returns true if the nuclide has coherrent elastic scattering.
   */
  bool has_coherent_elastic() const { return has_coherent_elastic_; }

  /**
   * @brief Returns true if the nuclide has incoherrent elastic scattering.
   */
  bool has_incoherent_elastic() const { return has_incoherent_elastic_; }

  /**
   * @brief Returns a STTSLReaction reference to the coherent elastic scattering
   *        data.
   */
  const STTSLReaction& coherent_elastic() const { return *coherent_elastic_; }

  /**
   * @brief Returns a STTSLReaction reference to the incoherent elastic
   *        scattering data.
   */
  const STTSLReaction& incoherent_elastic() const {
    return *incoherent_elastic_;
  }

  /**
   * @brief Returns a STTSLReaction reference to the incoherent inelastic data
   *        scattering data.
   */
  const STTSLReaction& incoherent_inelastic() const {
    return *incoherent_inelastic_;
  }

  /**
   * @breif Returns the total thermal scattering cross section.
   * @param E Energy at which to evaluate the cross section.
   */
  double xs(double E) const {
    double ii = incoherent_inelastic_->xs(E);
    double ie = incoherent_elastic_->xs(E);
    double ce = coherent_elastic_->xs(E);
    return ii + ie + ce;
  }

 private:
  ZAID zaid_;
  double awr_;
  double temperature_;
  bool has_coherent_elastic_;
  bool has_incoherent_elastic_;

  std::shared_ptr<STTSLReaction> coherent_elastic_;
  std::shared_ptr<STTSLReaction> incoherent_elastic_;
  std::shared_ptr<STIncoherentInelastic> incoherent_inelastic_;
};
}  // namespace pndl

#endif
