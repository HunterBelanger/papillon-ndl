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
#ifndef PANGLOS_COHERENT_ELASTIC_H
#define PANGLOS_COHERENT_ELASTIC_H

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
 * @brief This class contains all information for Coherent Elastic scattering.
 *        It is able to calculate the cross section at any temperature.
 */
class CoherentElastic {
 public:
  /**
   * @param ce Coherent Elastic information from ENDFtk.
   */
  CoherentElastic(const section::Type<7, 2>::CoherentElastic& ce);

  /**
   * @brief Returns the cross section for coherent elastic scattering.
   *        This therefore returns
            \f[
                \sigma_{CE}(T, E_{\text{in}}) =
                \frac{1}{E_{\text{in}}}
                \sum_{i}^{E_i < E_{\text{in}}}
                s_i(T)
            \f]

   * @param T Temperature in Kelvin.
   * @param Ein Incident energy in eV.
   */
  double xs(double T, double Ein) const;

  /**
   * @brief Returns a reference to a vector containing all of the Bragg edges
   *        with energy units of eV.
   */
  const std::vector<double>& bragg_edges() const { return bragg_edges_; }

  /**
   * @brief Returns a reference to a vector contianing the cummulative sum of
   *        all of the structure factors, for each Bragg edge, at a given
   *        temperature index.
   * @param i Indext of the desired temperature.
   */
  const std::vector<double>& structure_factors(std::size_t i) const {
    return structure_factor_sums_[i];
  }

  /**
   * @brief Interpolates temperature dependent structure factors to the desired
   *        temperature.
   * @param T Temperature for cummulative sum of structure factors in K.
   */
  std::vector<double> interpolate_structure_factors(double T) const;

  /**
   * @breif Returns true if the structure factors are temperature dependent.
   */
  bool temperature_dependent() const { return temperatures_.size() > 1; }

  /**
   * @brief Returns a reference to the vector of all tabulated temperatures.
   */
  const std::vector<double>& temperatures() { return temperatures_; }

  /**
   * @brief Returns a reference to the vector of the Interpolator instances,
   *        to interpolate between different temperatures. If only one
   *        temperature is provided, this vector is empty.
   */
  const std::vector<Interpolator>& temperature_interpolators() {
    return temp_interps_;
  }

 private:
  std::vector<double> bragg_edges_;
  std::vector<std::vector<double>> structure_factor_sums_;
  std::vector<double> temperatures_;
  std::vector<Interpolator> temp_interps_;
};

#endif
