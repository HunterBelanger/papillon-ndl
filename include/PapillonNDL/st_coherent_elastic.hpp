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
#ifndef PAPILLON_NDL_ST_COHERENT_ELASTIC_H
#define PAPILLON_NDL_ST_COHERENT_ELASTIC_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <algorithm>

namespace pndl {

/**
 * @brief Holds the Coherent Elastic scattering data for a single nuclide
 *        at a single temperature.
 */
class STCoherentElastic {
 public:
  /**
   * @param ace ACE file which contains thermal scattering law.
   */
  STCoherentElastic(const ACE& ace);
  ~STCoherentElastic() = default;

  /**
   * @brief Evaluates the incoherent inelastic scattering cross section
   *        at energy E.
   * @param E Incident energy at which to evaluate the cross section in MeV.
   */
  double xs(double E) const {
    if (E > bragg_edges_.front() && E < bragg_edges_.back()) {
      // Get index for lower bragg edge
      auto Eit = std::lower_bound(bragg_edges_.begin(), bragg_edges_.end(), E);
      size_t l = std::distance(bragg_edges_.begin(), Eit) - 1;
      return structure_factor_sum_[l] / E;
    } else if (E < bragg_edges_.front()) {
      return 0.;
    } else {
      return structure_factor_sum_.back() / E;
    }
  }

  /**
   * @brief Sample the angle-energy distribution.
   * @param E_in Incident energy in MeV.
   */
  AngleEnergyPacket sample_angle_energy(double E_in) const {
    // Get Bragg edge of scatter
    double Ei = 0.;
    if (E_in > bragg_edges_.front() && E_in < bragg_edges_.back()) {
      // Get index for lower bragg edge
      auto Eit =
          std::lower_bound(bragg_edges_.begin(), bragg_edges_.end(), E_in);
      Eit--;
      Ei = *Eit;
    } else if (E_in < bragg_edges_.front()) {
      Ei = 0.;
    } else {
      Ei = bragg_edges_.back();
    }

    double mu = 1. - (2. * Ei / E_in);

    return {mu, E_in};
  }

  /**
   * @brief Returns the vector of Bragg edges.
   */
  const std::vector<double>& bragg_edges() const { return bragg_edges_; }

  /**
   * @brief Returns the vector of the sum of structure factors.
   */
  const std::vector<double>& structure_factor_sum() const {
    return structure_factor_sum_;
  }

 private:
  std::vector<double> bragg_edges_;
  std::vector<double> structure_factor_sum_;
};

}  // namespace pndl

#endif