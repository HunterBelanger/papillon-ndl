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
#ifndef PAPILLON_NDL_CM_DISTRIBUTION_H
#define PAPILLON_NDL_CM_DISTRIBUTION_H

#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/frame.hpp>
#include <memory>
#include <vector>

/**
 * @file
 * @author Hunter Belanger
 */

namespace pndl {

/**
 * @brief A dsitribution which for which the data is provided in the center of
 *        mass frame.
 */
class CMDistribution : public AngleEnergy {
 public:
  /**
   * @param A Atomic weight ratio of the nuclide.
   * @param Q The Q-value of the reaction.
   * @param distribution Pointer to the distribution object in the center of
   * mass frame.
   *
   */
  CMDistribution(double A, double Q, std::shared_ptr<AngleEnergy> distribution)
      : awr_(A), q_(Q), distribution_(distribution) {}

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  std::optional<double> angle_pdf(double E_in, double mu) const override final;

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override final;
  /**
   * @brief Returns the distribution in the Center of Mass frame.
   */
  const AngleEnergy& distribution() const { return *distribution_; }

  /**
   * @brief Returns the nuclide Atomic Weight Ratio.
   */
  double awr() const { return awr_; }

  /**
   * @brief Returns the Q-value of the reaction.
   */
  double q() const { return q_; }

 private:
  double awr_, q_;
  std::shared_ptr<AngleEnergy> distribution_;
};

}  // namespace pndl

#endif
