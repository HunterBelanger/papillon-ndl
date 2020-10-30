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
#ifndef PAPILLON_NDL_CROSS_SECTION_H
#define PAPILLON_NDL_CROSS_SECTION_H

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_grid.hpp>
#include <PapillonNDL/shared_span.hpp>

namespace pndl {

class CrossSection {
 public:
  CrossSection() : energy_values_({0.}), values_(), index_(0) {}
  CrossSection(const ACE& ace, size_t i, const EnergyGrid& E_grid,
               bool get_index = true);
  ~CrossSection() = default;

  double operator[](size_t i) const {
    if (i < index_)
      return values_.front();
    else if (i >= index_ + values_.size())
      return values_.back();

    return values_[i - index_];
  }

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

  uint32_t index() const;

  size_t size() const;
  double value(size_t i) const;
  double energy(size_t i) const;

  const std::vector<float>& values() const;
  std::vector<float> energies() const;

 private:
  shared_span<float> energy_values_;
  std::vector<float> values_;
  uint32_t index_;
};

}  // namespace pndl

#endif
