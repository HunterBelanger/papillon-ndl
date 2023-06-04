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
#ifndef PAPILLON_NDL_ST_INCOHERENT_ELASTIC_H
#define PAPILLON_NDL_ST_INCOHERENT_ELASTIC_H

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
 *        at a single temperature, according to the custom Panglos ACE format.
 */
class STIncoherentElastic : public STTSLReaction {
 public:
  /**
   * @param ace ACE file which contains thermal scattering law.
   */
  STIncoherentElastic(const ACE& ace);
  ~STIncoherentElastic() = default;

  double xs(double E) const override final {
    if (xs_ < 0.) {
      return 0.;
    }

    return 0.5 * xs_ * ((1. - std::exp(-4. * E * W_)) / (2. * E * W_));
  }

  AngleEnergyPacket sample_angle_energy(
      double E_in, const std::function<double()>& rng) const override final {
    if (xs_ < 0.) {
      std::string mssg =
          "Incoherent elastic scattering is not possible. Cannot sample "
          "distribution.";
      throw PNDLException(mssg);
    }

    const double xi = rng();
    const double EW2 = 2. * E_in * W_;
    double mu = (std::log(xi * (std::exp(2. * EW2) - 1.) + 1.) / EW2) - 1.;

    if (std::abs(mu) > 1.) {
      if (mu > 1.)
        mu = 1.;
      else if (mu < -1.)
        mu = -1.;
    }

    return {mu, E_in};
  }

  std::optional<double> angle_pdf(double E_in, double mu) const override final {
    const double EW2 = 2. * E_in * W_;
    const double exp_2EW = std::exp(EW2);
    const double exp_neg2EW = 1. / exp_2EW;
    const double c = (1. / EW2) * (exp_2EW - exp_neg2EW);

    return c * std::exp(EW2 * mu);
  }

  std::optional<double> pdf(double /*E_in*/, double /*mu*/,
                            double /*E_out*/) const override final {
    return std::nullopt;
  }

  /**
   * @brief Returns the charicteristic bound cross section.
   */
  double bound_xs() const { return xs_; }

  /**
   * @brief Returns the Debye-Waller integral divided by the atomic mass.
   */
  double W() const { return W_; }

 private:
  double xs_;
  double W_;
};

}  // namespace pndl

#endif
