/*
 * Copyright 2020, Hunter Belanger
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

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/shared_span.hpp>
#include <cmath>
#include <vector>

namespace pndl {

class EnergyGrid {
 public:
  EnergyGrid(const ACE& ace, uint32_t NBINS = 8192);
  ~EnergyGrid() = default;

  double operator[](size_t i) const;
  size_t size() const;
  const shared_span<float>& grid() const;
  double min_energy() const;
  double max_energy() const;

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
                                 energy_values_.begin() + hi_indx, E) - energy_values_.begin() - 1;
    
    return ind;
  }

 private:
  shared_span<float> energy_values_;
  std::vector<uint32_t> bin_pointers_;
  double u_min, du;
};

}  // namespace pndl

#endif
