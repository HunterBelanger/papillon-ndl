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
#ifndef PAPILLON_NDL_ENERGY_GRID_H
#define PAPILLON_NDL_ENERGY_GRID_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <cmath>
#include <vector>

namespace pndl {

/**
 * @brief Holds the hashed energy grid of a Nuclide.
 */
class EnergyGrid {
 public:
  /**
   * @param ace ACE file from which to take the energy grid.
   * @param NBINS Number of bins to hash the energy grid into. The
   *              default value is 8192, which is the number of bins
   *              used by MCNP.
   */
  EnergyGrid(const ACE& ace, uint32_t NBINS = 8192);

  /**
   * @param energy Vector of all points in energy grid (sorted).
   * @param NBINS Number of bins to hash the energy grid into. The
   *              default value is 8192, which is the number of bins
   *              used by MCNP.
   */
  EnergyGrid(const std::vector<double>& energy, uint32_t NBINS = 8192);

  ~EnergyGrid() = default;

  /**
   * @brief Returns the ith energy in the grid in MeV.
   * @param i Index into energy grid.
   */
  double operator[](size_t i) const { return energy_values_[i]; }

  /**
   * @brief Number of points in the complete energy grid.
   */
  size_t size() const { return energy_values_.size(); }

  /**
   * @brief Returns a reference to the energy grid.
   */
  const std::vector<double>& grid() const { return energy_values_; }

  /**
   * @brief Returns the lowest energy in the grid.
   */
  double min_energy() const { return energy_values_.front(); }

  /**
   * @brief Returns the highest energy in the grid.
   */
  double max_energy() const { return energy_values_.back(); }

  /**
   * @brief Finds the interpolation index for a given energy, using the
   *        hashing algorithm for speed.
   * @param E Energy for which to find the index.
   */
  size_t get_lower_index(double E) const {
    if (E <= energy_values_.front()) {
      return 0;
    } else if (E >= energy_values_.back()) {
      return energy_values_.size() - 1;
    }

    // Get current bin
    uint32_t bin = static_cast<uint32_t>((std::log(E) - u_min) / du);

    // lower search index
    uint32_t low_indx = bin_pointers_[bin];
    uint32_t hi_indx = bin_pointers_[bin + 1] + 1;

    size_t ind = std::lower_bound(energy_values_.begin() + low_indx,
                                  energy_values_.begin() + hi_indx, E) -
                 energy_values_.begin() - 1;

    return ind;
  }

  /**
   * @brief Re-hashes the energy grid to specified number of pointers.
   * @param NBINS Number of bins to hash the energy grid into.
   */
  void hash_energy_grid(uint32_t NBINS);

 private:
  std::vector<double> energy_values_;
  std::vector<uint32_t> bin_pointers_;
  double u_min, du;
};

}  // namespace pndl

#endif
