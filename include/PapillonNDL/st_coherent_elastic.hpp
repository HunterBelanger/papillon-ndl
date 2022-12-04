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
#ifndef PAPILLON_NDL_ST_COHERENT_ELASTIC_H
#define PAPILLON_NDL_ST_COHERENT_ELASTIC_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/st_tsl_reaction.hpp>
#include <algorithm>
#include <iterator>
#include <optional>

namespace pndl {

/**
 * @brief Holds the Coherent Elastic scattering data for a single nuclide
 *        at a single temperature.
 */
class STCoherentElastic : public STTSLReaction {
 public:
  /**
   * @param ace ACE file which contains thermal scattering law.
   */
  STCoherentElastic(const ACE& ace);
  ~STCoherentElastic() = default;

  double xs(double E) const override final {
    if (bragg_edges_.size() == 0) return 0.;

    if (E > bragg_edges_.front() && E < bragg_edges_.back()) {
      // Get index for lower bragg edge
      auto Eit = std::lower_bound(bragg_edges_.begin(), bragg_edges_.end(), E);
      std::size_t l = static_cast<std::size_t>(
          std::distance(bragg_edges_.begin(), Eit) - 1);
      return structure_factor_sum_[l] / E;
    } else if (E < bragg_edges_.front()) {
      return 0.;
    } else {
      return structure_factor_sum_.back() / E;
    }
  }

  AngleEnergyPacket sample_angle_energy(
      double E_in, const std::function<double()>& rng) const override final {
    if (bragg_edges_.size() == 0) {
      std::string mssg =
          "Coherent elastic scattering is not possible. Cannot sample "
          "distribution.";
      throw PNDLException(mssg);
    }

    if (E_in > bragg_edges_.front()) {
      // Get index for lower bragg edge
      auto Eit =
          std::lower_bound(bragg_edges_.begin(), bragg_edges_.end(), E_in);
      std::size_t l = static_cast<std::size_t>(
          std::distance(bragg_edges_.begin(), Eit) - 1);

      // Sample which Bragg edge off of which we will scatter.
      double Prob = rng() * structure_factor_sum_[l];
      auto Sit = std::lower_bound(
          structure_factor_sum_.begin(),
          structure_factor_sum_.begin() + static_cast<std::ptrdiff_t>(l), Prob);
      std::size_t Si = static_cast<std::size_t>(
          std::distance(structure_factor_sum_.begin(), Sit));
      double E_bragg = bragg_edges_[Si];

      // Calculate the cosine of the scattering angle.
      double mu = 1. - (2. * E_bragg / E_in);

      return {mu, E_in};
    } else {
      // When E_in <= E_0, the xs is 0, so we shouldn't actually be sampling
      // this distribution. We will indicated this through a forward scatter
      // with no change in energy.
      return {1., E_in};
    }
  }

  std::optional<double> angle_pdf(double /*E_in*/,
                                  double /*mu*/) const override final {
    return std::nullopt;
  }

  std::optional<double> pdf(double /*E_in*/, double /*mu*/,
                            double /*E_out*/) const override final {
    return std::nullopt;
  }

  /**
   * @brief Returns the vector of Bragg edges.
   */
  const std::vector<double>& bragg_edges() const { return bragg_edges_; }

  /**
   * @brief Returns the vector of the sum of structure factors.
   */
  const std::vector<double>& structure_factor_sum() const {
    return structure_factor_sum_;
  }

 private:
  std::vector<double> bragg_edges_;
  std::vector<double> structure_factor_sum_;
};

}  // namespace pndl

#endif
