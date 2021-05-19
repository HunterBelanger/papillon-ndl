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
#include <PapillonNDL/equiprobable_angle_bins.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>

namespace pndl {

EquiprobableAngleBins::EquiprobableAngleBins(const ACE& ace, std::size_t i)
    : bounds_(ace.xss(i, NBOUNDS)) {
  if (!std::is_sorted(bounds_.begin(), bounds_.end())) {
    std::string mssg =
        "EquiprobableAngleBins::EquiprobableAngleBins: Bin bounds are not "
        "sorted. Index of EquiprobableAngleBins in XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (bounds_[0] < -1.) {
    std::string mssg =
        "EquiprobableAngleBins::EquiprobableAngleBins: Lowest bin bound is "
        "less than -1. Index of EquiprobableAngleBins in XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (bounds_[NBOUNDS - 1] > 1.) {
    std::string mssg =
        "EquiprobableAngleBins::EquiprobableAngleBins: Highest bin bound is "
        "more than 1. Index of EquiprobableAngleBins in XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }
}

EquiprobableAngleBins::EquiprobableAngleBins(const std::vector<double>& bounds)
    : bounds_(bounds) {
  if (bounds_.size() != 33) {
    std::string mssg =
        "EquiprobableAngleBins::EquiprobableAngleBins: Must provide 33 bin "
        "boundaries.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (!std::is_sorted(bounds_.begin(), bounds_.end())) {
    std::string mssg =
        "EquiprobableAngleBins::EquiprobableAngleBins: Bin bounds are not "
        "sorted.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (bounds_[0] < -1.) {
    std::string mssg =
        "EquiprobableAngleBins::EquiprobableAngleBins: Lowest bin bound is "
        "less than -1.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (bounds_[NBOUNDS - 1] > 1.) {
    std::string mssg =
        "EquiprobableAngleBins::EquiprobableAngleBins: Highest bin bound is "
        "more than 1.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }
}

double EquiprobableAngleBins::sample_mu(double xi) const {
  std::size_t bin =
      static_cast<std::size_t>(std::floor(static_cast<double>(NBOUNDS) * xi));
  if (bin == NBOUNDS) bin--;
  double C_b = bin * P_BIN;
  double mu_low = bounds_[bin];
  double mu = ((xi - C_b) / P_BIN) + mu_low;

  if (std::abs(mu) > 1.) mu = std::copysign(1, mu);

  return mu;
}

double EquiprobableAngleBins::pdf(double mu) const {
  if (mu < bounds_.front() || mu > bounds_.back()) return 0.;

  std::size_t bin = 0;
  for (std::size_t i = 0; i < NBOUNDS - 1; i++) {
    if (bounds_[i] <= mu && bounds_[i + 1] >= mu) {
      bin = i;
      break;
    }
  }

  double mu_low = bounds_[bin];
  double mu_hi = bounds_[bin + 1];
  return P_BIN / (mu_hi - mu_low);
}

}  // namespace pndl
