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
#ifndef PANGLOS_INCOHERENT_INELASTIC_H
#define PANGLOS_INCOHERENT_INELASTIC_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <ENDFtk/file/7.hpp>
#include <ENDFtk/section/7.hpp>

#include "tabulated_sab.hpp"
using namespace njoy::ENDFtk;

#include <memory>
#include <vector>

/**
 * @brief This class contains all information for Incoherent Inelastic
 *        scattering. It is able to calculate the double differential and
 *        integral scattering cross sections, for any of the tabulated
 *        temperatures provided in the the evaluation.
 */
class IncoherentInelastic {
 public:
  /**
   * @param mt4 Information for MF7 MT4 from ENDFtk.
   */
  IncoherentInelastic(const section::Type<7, 4>& mt4);

  /**
   * @brief Returns the double differential cross section for incoherent
   *        inelastic scattering. This therefore returns
            \f[
                \sigma_{II}(E_{\text{in}}\rightarrow E_{\text{out}}, \mu) =
                \frac{A kT \sigma_b}{4 E_{\text{in}}}
   e^{-\beta/2}S(\alpha,\beta) \f]
   * @param Ti Temperature index.
   * @param Ein Incident energy in eV.
   * @param Eout Exit energy in eV.
   * @param mu Cosine of scattering angle, in interval [-1, 1].
   */
  double ddxs(std::size_t Ti, double Ein, double Eout, double mu) const;

  /**
   * @brief Returns the cross section for incoherent inelastic scattering.
   *        This therefore returns
            \f[
                \sigma_{II}(E_{\text{in}}) =
                \frac{A kT \sigma_b}{4 E_{\text{in}}}
                \int_{\beta_{\text{min}}}^{\beta_{\text{max}}}
                \int_{\alpha_{\text{min}}}^{\alpha_{\text{max}}}
                e^{-\beta/2}S(\alpha,\beta)
                d\alpha d\beta
            \f]

   * @param Ti Temperature index.
   * @param Ein Incident energy in eV.
   */
  double xs(std::size_t Ti, double Ein) const;

  /**
   * @brief Returns reference to vector of all temperatures for which
   *        \f$S(\alpha,\beta)\f$ is tabulated.
   */
  const std::vector<double>& temperatures() const { return sab_temps_; }

  /**
   * @brief Returns a reference to a \f$S(\alpha,\beta)\f$ law for the given
   *        temperature index.
   * @param Ti Temerature index.
   */
  const TabulatedSab& sab(std::size_t Ti) const { return *sab_[Ti]; }

  /**
   * @brief Returns the minimum energy for the thermal scattering law data.
   */
  double Emin() const { return Emin_; }

  /**
   * @brief Returns the maximum energy for the thermal scattering law data.
   */
  double Emax() const { return Emax_; }

  /**
   * @brief Returns the bound cross section for a single istope of the principal
   *        scatterer, in barns.
   */
  double bound_xs() const { return bound_xs_; }

  /**
   * @brief Returns the atomic weight ratio of the principal scatterer.
   */
  double awr() const { return awr_; }

 private:
  std::vector<std::unique_ptr<TabulatedSab>> sab_;
  std::vector<double> sab_temps_;
  double awr_;
  double bound_xs_;
  double Emin_;
  double Emax_;

  void setup_ii(const section::Type<7, 4>& mt4);
};

struct AlphaDistribution {
  std::vector<double> alpha, pdf, cdf;
};

struct BetaDistribution {
  std::vector<double> beta, pdf, cdf;
  std::vector<AlphaDistribution> alpha;
};

struct LinearizedIncoherentInelastic {
  std::vector<double> egrid, xs;
  std::vector<BetaDistribution> beta;
};

LinearizedIncoherentInelastic linearize_ii(const IncoherentInelastic& ii,
                                           std::size_t Ti,
                                           bool pedantic);

#endif
