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
#ifndef PAPILLON_NDL_REACTION_BASE_H
#define PAPILLON_NDL_REACTION_BASE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/frame.hpp>
#include <PapillonNDL/function_1d.hpp>
#include <memory>

namespace pndl {

/**
 * @brief Holds the product distributions for a single MT.
 */
class ReactionBase {
 public:
  ReactionBase(const ReactionBase&) = default;

  /**
   * @brief Returns the MT of the reaction.
   */
  uint32_t mt() const { return mt_; }

  /**
   * @brief Returns the Q-value of the reaction.
   */
  double q() const { return q_; }

  /**
   * @brief Returns the threshold energy for the reaction.
   */
  double threshold() const { return threshold_; }

  /**
   * @brief Returns the function for the reaction yield.
   */
  const Function1D& yield() const { return *yield_; }

  /**
   * @brief Samples and angle and energy from the neutron reaction
   *        product distribution.
   * @param E_in Incident energy in MeV.
   * @param rng Random number generation function.
   */
  AngleEnergyPacket sample_neutron_angle_energy(
      double E_in, std::function<double()> rng) const {
    if (E_in < threshold_) return {0., 0.};

    return neutron_distribution_->sample_angle_energy(E_in, rng);
  }

  /**
   * @brief Returns the distribution for neutron reaction products.
   */
  const AngleEnergy& neutron_distribution() const {
    return *neutron_distribution_;
  }

 protected:
  uint32_t mt_;
  double q_;
  double awr_;
  double threshold_;
  std::shared_ptr<Function1D> yield_;
  std::shared_ptr<AngleEnergy> neutron_distribution_;

  /**
   * @param ace ACE file to take reaction from.
   * @param indx Reaction index in the MT array.
   */
  ReactionBase(const ACE& ace, std::size_t indx);

  // Private helper methods
  void load_neutron_distributions(
      const ACE& ace, std::size_t indx,
      std::vector<std::shared_ptr<AngleEnergy>>& distributions,
      std::vector<std::shared_ptr<Tabulated1D>>& probabilities);
};

}  // namespace pndl

#endif
