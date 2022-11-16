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
#ifndef PANGLOS_INCOHERENT_ELASTIC_H
#define PANGLOS_INCOHERENT_ELASTIC_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <ENDFtk/file/7.hpp>
#include <ENDFtk/section/7.hpp>
using namespace njoy::ENDFtk;

#include <memory>
#include <vector>

#include "interpolator.hpp"

/**
 * @brief This class contains all information for Incoherent Elastic scattering.
 *        It is able to calculate the single differential and integral
 *        scattering cross sections, at any temperature.
 */
class IncoherentElastic {
 public:
  /**
   * @param ie Incoherent Elastic information from ENDFtk.
   */
  IncoherentElastic(const section::Type<7, 2>::IncoherentElastic& ie);

  /**
   * @brief Returns the single differential cross section for incoherent
   *        elastic scattering. This therefore returns
            \f[
                \sigma_{IE}(T, E_{\text{in}}, \mu) =
                \frac{\sigma_b}{4\pi}e^{-2E_{\text{in}}W'(T)(1-\mu)}
            \f]
   * @param T Temperature in Kevlin.
   * @param Ein Incident energy in eV.
   * @param mu Cosine of scattering angle, in interval [-1, 1].
   */
  double dxs(double T, double Ein, double mu) const;

  /**
   * @brief Returns the cross section for incoherent elastic scattering.
   *        This therefore returns
            \f[
                \sigma_{IE}(T, E_{\text{in}}) =
                \frac{\sigma_b}{2}
                \left(
                   \frac{1-e^{-4E_{\text{in}}W'(T)}}
                        {2E_{\text{in}}W'(T)}
                \right)
            \f]

   * @param T Temperature in Kelvin.
   * @param Ein Incident energy in eV.
   */
  double xs(double T, double Ein) const;

  /**
   * @brief Returns the bound cross section.
   */
  double bound_xs() const { return bound_xs_; }

  /**
   * @brief Returns the reference to the tabulated Debye-Waller integral divided
   *        by the atomic mass.
   */
  const Tab1& W() const { return W_; }

 private:
  Tab1 W_;
  double bound_xs_;
};

#endif
