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
#ifndef PAPILLON_NDL_PCTABLE_H
#define PAPILLON_NDL_PCTABLE_H

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/interpolation.hpp>
#include <cmath>
#include <vector>

namespace pndl {

class PCTable {
 public:
  PCTable(const ACE& ace, size_t i, double normalization = 1.);
  ~PCTable() = default;

  double sample_value(double xi) const {
    if (xi < 0. || xi > 1.) {
      throw PNDLException("PCTable: Invalid value for xi provided.", __FILE__, __LINE__);
    }

    auto cdf_it = std::lower_bound(cdf_.begin(), cdf_.end(), xi);
    size_t l = std::distance(cdf_.begin(), cdf_it);
    if (xi == *cdf_it) return values_[l];

    l--;

    if (interp_ == Interpolation::Histogram) return histogram_interp(xi, l);

    return linear_interp(xi, l);
  }

  double min_value() const { return values_.front(); }
  double max_value() const { return values_.back(); }

  size_t size() const { return values_.size(); }
  const std::vector<double>& values() const { return values_; }
  const std::vector<double>& pdf() const { return pdf_; }
  const std::vector<double>& cdf() const { return cdf_; }
  Interpolation interpolation() const { return interp_; }

 private:
  std::vector<double> values_;
  std::vector<double> pdf_;
  std::vector<double> cdf_;
  Interpolation interp_;

  double histogram_interp(double xi, size_t l) const {
    return values_[l] + ((xi - cdf_[l]) / pdf_[l]);
  }

  double linear_interp(double xi, size_t l) const {
    double m = (pdf_[l + 1] - pdf_[l]) / (values_[l + 1] - values_[l]);
    return values_[l] +
           (1. / m) * (std::sqrt(pdf_[l] * pdf_[l] + 2. * m * (xi - cdf_[l])) -
                       pdf_[l]);
  }
};

}  // namespace pndl

#endif
