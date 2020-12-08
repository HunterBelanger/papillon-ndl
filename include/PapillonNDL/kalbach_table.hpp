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
#ifndef PAPILLON_NDL_KALBACH_TABLE_H
#define PAPILLON_NDL_KALBACH_TABLE_H

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/interpolation.hpp>
#include <algorithm>
#include <cmath>

namespace pndl {

class KalbachTable {
 public:
  KalbachTable(const ACE& ace, size_t i);
  ~KalbachTable() = default;

  double sample_energy(double xi) const {
    if (xi < 0. || xi > 1.) {
      throw std::runtime_error("KalbachTable: Invalid value for xi provided");
    }

    auto cdf_it = std::lower_bound(cdf_.begin(), cdf_.end(), xi);
    size_t l = std::distance(cdf_.begin(), cdf_it);
    if (xi == *cdf_it) return energy_[l];
    l--;

    if (interp_ == Interpolation::Histogram)
      return histogram_interp_energy(xi, l);

    return linear_interp_energy(xi, l);
  }

  double min_energy() const { return energy_.front(); }

  double max_energy() const { return energy_.back(); }

  double R(double E) const {
    if (E <= energy_.front())
      return R_.front();
    else if (E >= energy_.back())
      return R_.back();
    else {
      auto E_it = std::lower_bound(energy_.begin(), energy_.end(), E);
      size_t l = std::distance(energy_.begin(), E_it) - 1;

      return interpolate(E, energy_[l], R_[l], energy_[l + 1], R_[l + 1],
                         interp_);
    }
  }

  double A(double E) const {
    if (E <= energy_.front())
      return A_.front();
    else if (E >= energy_.back())
      return A_.back();
    else {
      auto E_it = std::lower_bound(energy_.begin(), energy_.end(), E);
      size_t l = std::distance(energy_.begin(), E_it) - 1;

      return interpolate(E, energy_[l], A_[l], energy_[l + 1], A_[l + 1],
                         interp_);
    }
  }

  const std::vector<double>& energy() const { return energy_; }
  const std::vector<double>& pdf() const { return pdf_; }
  const std::vector<double>& cdf() const { return cdf_; }
  const std::vector<double>& R() const { return R_; }
  const std::vector<double>& A() const { return A_; }
  Interpolation interpolation() const { return interp_; }

 private:
  std::vector<double> energy_;
  std::vector<double> pdf_;
  std::vector<double> cdf_;
  std::vector<double> R_;
  std::vector<double> A_;
  Interpolation interp_;

  double histogram_interp_energy(double xi, size_t l) const {
    return energy_[l] + ((xi - cdf_[l]) / pdf_[l]);
  }

  double linear_interp_energy(double xi, size_t l) const {
    double m = (pdf_[l + 1] - pdf_[l]) / (energy_[l + 1] - energy_[l]);
    return energy_[l] +
           (1. / m) * (std::sqrt(pdf_[l] * pdf_[l] + 2. * m * (xi - cdf_[l])) -
                       pdf_[l]);
  }
};

}  // namespace pndl

#endif
