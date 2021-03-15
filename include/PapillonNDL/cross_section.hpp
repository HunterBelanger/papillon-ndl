/*
 * Copyright 2021, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
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
#include <PapillonNDL/shared_span.hpp>

namespace pndl {

/**
 * @brief Contains the linearly interpolable cross section data for
 *        a single MT.
 */
class CrossSection {
 public:
  /**
   * @param ace ACE file to take the data from.
   * @param i Index in the XSS block where the cross section starts.
   * @param E_grid Energy grid associated with the cross section values.
   * @param get_index Flag to indicate wether the cross section values begin
   *                  at i, or if the energy grid index is at i. Default
   *                  value is true.
   */
  CrossSection(const ACE& ace, size_t i, const EnergyGrid& E_grid,
               bool get_index = true);
  ~CrossSection() = default;

  /**
   * @brief Returns value of the cross section at index relative to
   *        the associated energy grid.
   * @param i Index from associated energy grid.
   */
  double operator[](size_t i) const {
    if (i < index_)
      return values_.front();
    else if (i >= index_ + values_.size())
      return values_.back();

    return values_[i - index_];
  }

  /**
   * @brief Evaluates the cross section at a given energy. Uses
   *        bisection search.
   * @param E Energy to evaluate the cross section at.
   */
  double operator()(double E) const {
    if (E <= energy_values_.front())
      return values_.front();
    else if (E >= energy_values_.back())
      return values_.back();

    auto E_it =
        std::lower_bound(energy_values_.begin(), energy_values_.end(), E);
    size_t indx = std::distance(energy_values_.begin(), E_it) - 1;

    double sig_low = values_[indx];
    double sig_hi = values_[indx + 1];
    double E_low = energy_values_[indx];
    double E_hi = energy_values_[indx + 1];

    return ((E - E_low) / (E_hi - E_low)) * (sig_hi - sig_low) + sig_low;
  }

  /**
   * @brief Evaluates the cross section at a given energy, with the
   *        grid point already provided.
   * @param E Energy to evaluate the cross section at.
   * @param i Index of the points for interpolation in the frame of
   *          the energy grid.
   */
  double operator()(double E, size_t i) const {
    if (E <= energy_values_.front())
      return values_.front();
    else if (E >= energy_values_.back())
      return values_.back();

    // Transform index from global grid to local grid
    i -= index_;

    double sig_low = values_[i];
    double sig_hi = values_[i + 1];
    double E_low = energy_values_[i];
    double E_hi = energy_values_[i + 1];

    return ((E - E_low) / (E_hi - E_low)) * (sig_hi - sig_low) + sig_low;
  }

  /**
   * @brief Evaluates the cross section at a given energy. Uses
   *        bisection search.
   * @param E Energy to evaluate the cross section at.
   */
  double evaluate(double E) const {
    return this->operator()(E);
  }
  
  /**
   * @brief Evaluates the cross section at a given energy, with the
   *        grid point already provided.
   * @param E Energy to evaluate the cross section at.
   * @param i Index of the points for interpolation in the frame of
   *          the energy grid.
   */
  double evaluate(double E, size_t i) const {
    return this->operator()(E, i);
  }

  /**
   * @brief Returns index in the energy grid at which the cross section
   *        values begin.
   */
  uint32_t index() const;

  /**
   * @brief Number of points in the cross section.
   */
  size_t size() const;

  /**
   * @brief Returns the ith cross section value.
   */
  double xs(size_t i) const;

  /**
   * @brief Returns the ith energy value, which corresponds with
   *        the ith cross section value.
   */
  double energy(size_t i) const;

  /**
   * @brief Returns the cross section values as a vector of floats.
   */
  const std::vector<float>& xs() const;

  /**
   * @brief Returns a copy of the energy grid points for the cross section
   *        as a vector of floats.
   */
  std::vector<float> energy() const;

 private:
  shared_span<float> energy_values_;
  std::vector<float> values_;
  uint32_t index_;
};

}  // namespace pndl

#endif
