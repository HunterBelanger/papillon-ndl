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
#ifndef PAPILLON_NDL_ANGLE_TABLE_H
#define PAPILLON_NDL_ANGLE_TABLE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/angle_law.hpp>
#include <PapillonNDL/legendre.hpp>
#include <PapillonNDL/pctable.hpp>
#include <functional>

namespace pndl {

/**
 * @brief Angular distribution which is provided as tabulated PDF and CDF.
 */
class AngleTable : public AngleLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   */
  AngleTable(const ACE& ace, std::size_t i);

  /**
   * @param cosines Conines of scattering angle which are tabulated.
   * @param pdf The Probability Density Function for the provided values.
   * @param cdf The Cumulative Density Function for the provided values.
   * @param interp Interpolation rule for the data. May be either
   *               Histogram or LinLin.
   */
  AngleTable(const std::vector<double>& cosines, const std::vector<double>& pdf,
             const std::vector<double>& cdf, Interpolation interp);

  /**
   * @param legendre Legendre distribution which will be linearized to create an
   *                 AngleTable.
   */
  AngleTable(const Legendre& legendre);

  /**
   * @param table PCTable contianing the PDF and CDF for the cosine
   *              distribution.
   */
  AngleTable(const PCTable& table);
  ~AngleTable() = default;

  double sample_mu(std::function<double()> rng) const override final;

  double pdf(double mu) const override final { return distribution_.pdf(mu); }

  /**
   * @brief Returns the number of points in the tabulated data.
   */
  std::size_t size() const { return distribution_.size(); }

  /**
   * @brief Returns the vector of the cosine points.
   */
  const std::vector<double>& cosines() const { return distribution_.values(); }

  /**
   * @brief Returns the vector of the PDF values.
   */
  const std::vector<double>& pdf() const { return distribution_.pdf(); }

  /**
   * @brief Returns the vector of the CDF values.
   */
  const std::vector<double>& cdf() const { return distribution_.cdf(); }

  /**
   * @brief Returns the type of interpolation used on the table
   *        (Histogram or LinLin).
   */
  Interpolation interpolation() const { return distribution_.interpolation(); }

 private:
  PCTable distribution_;
};

}  // namespace pndl

#endif
