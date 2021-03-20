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

EquiprobableAngleBins::EquiprobableAngleBins(const ACE& ace, size_t i)
    : bounds_(ace.xss(i, NBOUNDS)) {
  if (!std::is_sorted(bounds_.begin(), bounds_.end())) {
    std::string mssg =
        "EquiprobableAngleBins::EquiprobableAngleBins: Bin bounds are not "
        "sorted.\n";
    mssg += "Index of EquiprobableAngleBins in XSS block is " +
            std::to_string(i) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (bounds_[0] < -1.) {
    std::string mssg =
        "EquiprobableAngleBins::EquiprobableAngleBins: Lowest bin bound is "
        "less than -1.\n";
    mssg += "Index of EquiprobableAngleBins in XSS block is " +
            std::to_string(i) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (bounds_[NBOUNDS - 1] > 1.) {
    std::string mssg =
        "EquiprobableAngleBins::EquiprobableAngleBins: Highest bin bound is "
        "more than 1.\n";
    mssg += "Index of EquiprobableAngleBins in XSS block is " +
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
  size_t bin =
      static_cast<size_t>(std::floor(static_cast<double>(NBOUNDS) * xi));
  double C_b = bin * P_BIN;
  double mu_low = bounds_[bin];
  return ((xi - C_b) / P_BIN) + mu_low;
}

size_t EquiprobableAngleBins::size() const { return NBOUNDS; }

const std::vector<double>& EquiprobableAngleBins::bin_bounds() const {
  return bounds_;
}

}  // namespace pndl
