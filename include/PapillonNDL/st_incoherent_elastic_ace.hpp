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
#ifndef PAPILLON_NDL_ST_INCOHERENT_ELASTIC_ACE_H
#define PAPILLON_NDL_ST_INCOHERENT_ELASTIC_ACE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/st_tsl_reaction.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <algorithm>

namespace pndl {

/**
 * @brief Holds the Incoherent Elastic scattering data for a single nuclide
 *        at a single temperature, according to the ACE format.
 */
class STIncoherentElasticACE : public STTSLReaction {
 public:
  /**
   * @param ace ACE file which contains thermal scattering law.
   */
  STIncoherentElasticACE(const ACE& ace);
  ~STIncoherentElasticACE() = default;

  double xs(double E) const override final { return xs_->evaluate(E); }

  AngleEnergyPacket sample_angle_energy(
      double E_in, const std::function<double()>& rng) const override final {
    if (incoming_energy_.size() == 0) {
      std::string mssg =
          "Incoherent elastic scattering is not possible. Cannot sample "
          "distribution.";
      throw PNDLException(mssg);
    }

    // Get energy index
    auto Eit = std::lower_bound(incoming_energy_.begin(),
                                incoming_energy_.end(), E_in);
    size_t i = 0;
    double f = 0.;
    if (Eit == incoming_energy_.begin()) {
      i = 0;
      f = 0.;
    } else if (Eit == incoming_energy_.end()) {
      i = incoming_energy_.size() - 2;
      f = 1.;
    } else {
      i = std::distance(incoming_energy_.begin(), Eit) - 1;
      f = (E_in - incoming_energy_[i]) /
          (incoming_energy_[i + 1] - incoming_energy_[i]);
    }

    // Sample random index for cosine
    uint32_t j = static_cast<uint32_t>(Nmu * rng());

    double mu_prime =
        cosines_[i][j] + f * (cosines_[i + 1][j] - cosines_[i][j]);

    double mu_left = -1. - (mu_prime + 1.);
    if (j != 0) {
      mu_left = cosines_[i][j - 1] +
                f * (cosines_[i + 1][j - 1] - cosines_[i][j - 1]);
    }

    double mu_right = 1. - (mu_prime - 1.);
    if (j != Nmu - 1) {
      mu_right = cosines_[i][j + 1] +
                 f * (cosines_[i + 1][j + 1] - cosines_[i][j + 1]);
    }

    double mu = mu_prime + std::min(mu_prime - mu_left, mu_right - mu_prime) *
                               (rng() - 0.5);

    return {mu, E_in};
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
   * @brief Returns the cross section function.
   */
  const Tabulated1D& xs() const { return *xs_; }

  /**
   * @brief Returns vector to the incoming energy grid.
   */
  const std::vector<double>& incoming_energy() const {
    return incoming_energy_;
  }

  /**
   * @brief Returns array of discrete scattering cosines.
   */
  const std::vector<std::vector<double>>& cosines() const { return cosines_; }

 private:
  std::shared_ptr<Tabulated1D> xs_;
  uint32_t Nmu;
  std::vector<double> incoming_energy_;
  std::vector<std::vector<double>> cosines_;
};

}  // namespace pndl

#endif
