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
#include <PapillonNDL/energy_grid.hpp>
#include <algorithm>
#include <cmath>

#include "constants.hpp"

namespace pndl {

EnergyGrid::EnergyGrid(const ACE& ace, uint32_t NBINS)
    : energy_values_({0.}), bin_pointers_(), u_min(), du() {
  std::vector<float> raw_egrid = ace.xss<float>(ace.ESZ(), ace.nxs(2));
  for (auto& E : raw_egrid) E *= MEV_TO_EV;
  energy_values_ = shared_span<float>(raw_egrid.begin(), raw_egrid.end());

  // Generate pointers for lethargy bins
  u_min = std::log(energy_values_.front());
  double u_max = std::log(energy_values_.back());
  du = static_cast<double>(NBINS) / (u_max - u_min);

  bin_pointers_.reserve(NBINS);

  // Start by storing index to u_min which is 0
  bin_pointers_.push_back(0);

  // Get energy index for each lethargy bin bound
  size_t i = 0;
  for (size_t b = 1; b < NBINS; b++) {
    double u = u_min + b * du;
    double E = std::exp(u);

    auto it =
        std::lower_bound(energy_values_.begin() + i, energy_values_.end(), E);
    i = std::distance(energy_values_.begin(), it) - 1;

    bin_pointers_.push_back(i);
  }

  // Save pointer for last energy
  bin_pointers_.push_back(energy_values_.size() - 1);
}

double EnergyGrid::operator[](size_t i) const { return energy_values_[i]; }

size_t EnergyGrid::size() const { return energy_values_.size(); }

size_t EnergyGrid::get_lower_index(double E) {
  if (E <= energy_values_.front()) {
    return 0;
  } else if (E >= energy_values_.back()) {
    return energy_values_.size() - 1;
  }

  // Get current bin
  uint32_t bin = std::floor((std::log(E) - u_min) / du);

  // lower search index
  uint32_t low_indx = bin_pointers_[bin];
  uint32_t hi_indx = bin_pointers_[bin + 1] + 1;
  if (hi_indx > energy_values_.size() - 1) hi_indx = energy_values_.size() - 1;

  auto E_it = std::lower_bound(energy_values_.begin() + low_indx,
                               energy_values_.begin() + hi_indx, E);
  size_t ind = std::distance(energy_values_.begin(), E_it);
  if (*E_it != E) ind--;

  return ind;
}

const shared_span<float>& EnergyGrid::grid() const { return energy_values_; }

double EnergyGrid::min_energy() const { return energy_values_.front(); }

double EnergyGrid::max_energy() const { return energy_values_.back(); }

}  // namespace pndl
