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
#ifndef PAPILLON_NDL_CROSS_SECTION_H
#define PAPILLON_NDL_CROSS_SECTION_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_grid.hpp>
#include <memory>

namespace pndl {

/**
 * @brief Contains the linearly interpolable cross section data for
 *        a single MT.
 */
class CrossSection {
  // CrossSection instances should mainly be shared pointers in this
  // library. This is because sometimes certain cross sections aren't
  // present (such as photon production), and this is represented as a
  // nullptr. It also makes it easier to load a nuclide and just grab
  // the elastic scattering xs should someone want to do DBRC or
  // something like that.

 public:
  /**
   * @param ace ACE file to take the data from.
   * @param i Index in the XSS block where the cross section starts.
   * @param E_grid Pointer to EnergyGrid associated with the cross section
   * values.
   * @param get_index Flag to indicate wether the cross section values begin
   *                  at i, or if the energy grid index is at i. Default
   *                  value is true.
   */
  CrossSection(const ACE& ace, std::size_t i,
               std::shared_ptr<EnergyGrid> E_grid, bool get_index = true);

  /**
   * @param xs Vector containing the cross section values.
   * @param E_grid Pointer to EnergyGrid to use for the cross section.
   * @param index Starting index in the energy grid.
   */
  CrossSection(const std::vector<double>& xs,
               std::shared_ptr<EnergyGrid> E_grid, std::size_t index);

  /**
   * @param xs Value for the cross section at all points in the provided
   *           energy grid.
   * @param E_grid Pointer to EnergyGrid to use for the cross section.
   */
  CrossSection(double xs, std::shared_ptr<EnergyGrid> E_grid);

  ~CrossSection() = default;

  /**
   * @brief Returns value of the cross section at index relative to
   *        the associated energy grid.
   * @param i Index from associated energy grid.
   */
  double operator[](std::size_t i) const {
    if (single_value_) return values_->front();

    if (i < index_)
      return 0.;
    else if (i >= index_ + values_->size())
      return values_->back();

    return (*values_)[i - index_];
  }

  /**
   * @brief Evaluates the cross section at a given energy. Uses
   *        bisection search.
   * @param E Energy to evaluate the cross section at.
   */
  double operator()(double E) const {
    if (single_value_) {
      if (E < energy_grid_->min_energy()) return 0.;
      return values_->front();
    }

    if (E <= (*energy_grid_)[index_])
      return 0.;
    else if (E >= energy_grid_->max_energy())
      return values_->back();

    const auto& egrid = energy_grid_->grid();
    const auto erange_begin = egrid.begin() + index_;
    const auto erange_end = egrid.end();

    auto E_it = std::lower_bound(erange_begin, erange_end, E);
    std::size_t i = std::distance(erange_begin, E_it) - 1;

    double E_low = egrid[index_ + i];
    double E_hi = egrid[index_ + i + 1];
    double sig_low = (*values_)[i];
    double sig_hi = (*values_)[i + 1];

    return ((E - E_low) / (E_hi - E_low)) * (sig_hi - sig_low) + sig_low;
  }

  /**
   * @brief Evaluates the cross section at a given energy, with the
   *        grid point already provided.
   * @param E Energy to evaluate the cross section at.
   * @param i Index of the points for interpolation in the frame of
   *          the energy grid.
   */
  double operator()(double E, std::size_t i) const {
    if (single_value_) {
      if (E < energy_grid_->min_energy()) return 0.;
      return values_->front();
    }

    if (i < index_)
      return 0.;
    else if (i >= index_ + values_->size() - 1)
      return values_->back();

    // Transform index from global grid to local grid
    i -= index_;

    double E_low = (*energy_grid_)[index_ + i];
    double E_hi = (*energy_grid_)[index_ + i + 1];
    double sig_low = (*values_)[i];
    double sig_hi = (*values_)[i + 1];

    return ((E - E_low) / (E_hi - E_low)) * (sig_hi - sig_low) + sig_low;
  }

  /**
   * @brief Evaluates the cross section at a given energy, with the
   *        grid point already provided.
   * @param E Energy to evaluate the cross section at.
   * @param i Index of the points for interpolation in the frame of
   *          the energy grid.
   * @param El The value of the energy grid at index i.
   * @param Eh The value of the energy grid at index i+1.
   */
  double operator()(double E, std::size_t i, double El, double Eh) const {
    if (single_value_) {
      if (E < energy_grid_->min_energy()) return 0.;
      return values_->front();
    }

    if (i < index_)
      return 0.;
    else if (i >= index_ + values_->size() - 1)
      return values_->back();

    // Transform index from global grid to local grid
    i -= index_;

    double sig_low = (*values_)[i];
    double sig_hi = (*values_)[i + 1];

    return ((E - El) / (Eh - El)) * (sig_hi - sig_low) + sig_low;
  }

  /**
   * @brief Evaluates the cross section at a given energy. Uses
   *        bisection search.
   * @param E Energy to evaluate the cross section at.
   */
  double evaluate(double E) const { return this->operator()(E); }

  /**
   * @brief Evaluates the cross section at a given energy, with the
   *        grid point already provided.
   * @param E Energy to evaluate the cross section at.
   * @param i Index of the points for interpolation in the frame of
   *          the energy grid.
   */
  double evaluate(double E, std::size_t i) const {
    return this->operator()(E, i);
  }

  /**
   * @brief Evaluates the cross section at a given energy, with the
   *        grid point already provided.
   * @param E Energy to evaluate the cross section at.
   * @param i Index of the points for interpolation in the frame of
   *          the energy grid.
   * @param El The value of the energy grid at index i.
   * @param Eh The value of the energy grid at index i+1.
   */
  double evaluate(double E, std::size_t i, double El, double Eh) const {
    return this->operator()(E, i, El, Eh);
  }

  /**
   * @brief Returns index in the energy grid at which the cross section
   *        values begin.
   */
  std::size_t index() const { return index_; }

  /**
   * @brief Number of points in the cross section.
   */
  std::size_t size() const { return values_->size(); }

  /**
   * @brief Returns the ith cross section value.
   */
  double xs(std::size_t i) const { return (*values_)[i]; }

  /**
   * @brief Returns the ith energy value, which corresponds with
   *        the ith cross section value.
   */
  double energy(std::size_t i) const { return (*energy_grid_)[index_ + i]; }

  /**
   * @brief Returns the cross section values as a vector of floats.
   */
  const std::vector<double>& xs() const { return *values_; }

  /**
   * @brief Returns a copy of the energy grid points for the cross section
   *        as a vector of floats.
   */
  std::vector<double> energy() const;

 private:
  std::shared_ptr<EnergyGrid> energy_grid_;
  std::shared_ptr<std::vector<double>> values_;
  std::size_t index_;
  bool single_value_;
};

}  // namespace pndl

#endif
