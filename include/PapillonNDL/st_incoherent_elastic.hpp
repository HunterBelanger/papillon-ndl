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
#ifndef PAPILLON_NDL_ST_INCOHERENT_ELASTIC_H
#define PAPILLON_NDL_ST_INCOHERENT_ELASTIC_H

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
 * @brief Holds the Incoherent Elastic scattering data for a single nuclide
 *        at a single temperature.
 */
class STIncoherentElastic {
 public:
  /**
   * @param ace ACE file which contains thermal scattering law.
   */
  STIncoherentElastic(const ACE& ace);
  ~STIncoherentElastic() = default;

  /**
   * @brief Returns a pointer to the cross section function.
   */
  std::shared_ptr<Region1D> cross_section() const { return xs_; }

  /**
   * @brief Evaluates the incoherent inelastic scattering cross section
   *        at energy E.
   * @param E Incident energy at which to evaluate the cross section in MeV.
   */
  double xs(double E) const { return (*this->xs_)(E); }

  /**
   * @brief Sample the angle-energy distribution.
   * @param E_in Incident energy in MeV.
   */
  AngleEnergyPacket sample_angle_energy(double E_in,
                                        std::function<double()> rng) const {
    // Get energy index
    auto Eit = std::lower_bound(incoming_energy_.begin(),
                                incoming_energy_.end(), E_in);
    size_t i = 0;
    double f = 0.;
    if (Eit == incoming_energy_.begin()) {
      i = 0;
      f = 0.;
    } else if (Eit == incoming_energy_.end()) {
      i = incoming_energy_.size() - 2;
      f = 1.;
    } else {
      i = std::distance(incoming_energy_.begin(), Eit) - 1;
      f = (E_in - incoming_energy_[i]) /
          (incoming_energy_[i + 1] - incoming_energy_[i]);
    }

    // Sample random index for cosine
    uint32_t j = Nmu * rng();

    double mu_prime =
        cosines_[i][j] + f * (cosines_[i + 1][j] - cosines_[i][j]);

    double mu_left = -1. - (mu_prime + 1.);
    if (j != 0) {
      mu_left = cosines_[i][j - 1] +
                f * (cosines_[i + 1][j - 1] - cosines_[i][j - 1]);
    }

    double mu_right = 1. - (mu_prime - 1.);
    if (j != Nmu - 1) {
      mu_right = cosines_[i][j + 1] +
                 f * (cosines_[i + 1][j + 1] - cosines_[i][j + 1]);
    }

    double mu = mu_prime + std::min(mu_prime - mu_left, mu_prime + mu_right) *
                               (rng() - 0.5);

    return {mu, E_in};
  }

  /**
   * @brief Returns vector to the incoming energy grid.
   */
  const std::vector<double>& incoming_energy() const {
    return incoming_energy_;
  }

  /**
   * @brief Returns array of discrete scattering cosines.
   */
  const std::vector<std::vector<double>>& cosines() const { return cosines_; }

 private:
  std::shared_ptr<Region1D> xs_;
  uint32_t Nmu;
  std::vector<double> incoming_energy_;
  std::vector<std::vector<double>> cosines_;
};

}  // namespace pndl

#endif