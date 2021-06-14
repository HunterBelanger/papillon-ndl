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
#ifndef PAPILLON_NDL_MULTIPLE_DISTRIBUTION_H
#define PAPILLON_NDL_MULTIPLE_DISTRIBUTION_H

#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <vector>

/**
 * @file
 * @author Hunter Belanger
 */

namespace pndl {

/**
 * @brief A dsitribution which is composed of mutliple possible
 *        distributions, each with a tabulated probability.
 */
class MultipleDistribution : public AngleEnergy {
 public:
  MultipleDistribution(
      const std::vector<std::shared_ptr<AngleEnergy>>& distributions,
      const std::vector<std::shared_ptr<Tabulated1D>>& probabilities);

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  std::optional<double> angle_pdf(double E_in, double mu) const override final;

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override final;

  /**
   * @brief Returns the number of distributions for the reaction.
   */
  std::size_t size() const { return distributions_.size(); }

  /**
   * @brief Returns the ith distribution for the reaction.
   * @param i Index of distribution to fetch.
   */
  const AngleEnergy& distribution(std::size_t i) const {
    return *distributions_[i];
  }

  /**
   * @brief Returns the ith distribution's probability function.
   * @param i Index of distribution to fetch.
   */
  const Tabulated1D& probability(std::size_t i) const {
    return *probabilities_[i];
  }

 private:
  std::vector<std::shared_ptr<AngleEnergy>> distributions_;
  std::vector<std::shared_ptr<Tabulated1D>> probabilities_;
};

}  // namespace pndl

#endif
