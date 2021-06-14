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
#ifndef PAPILLON_NDL_DISCRETE_COSINES_ENERGIES_H
#define PAPILLON_NDL_DISCRETE_COSINES_ENERGIES_H

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>

/**
 * @file
 * @author Hunter Belanger
 */

namespace pndl {

/**
 * @brief Class which represents equiprobably and skewed discrete energy and
 *        discrete cosine distributions for incoherent inelastic scattering.
 */
class DiscreteCosinesEnergies : public AngleEnergy {
 public:
  /**
   * @param ace ACE file which contains thermal scattering law.
   */
  DiscreteCosinesEnergies(const ACE& ace);

  /**
   * @brief Struct to contain a discrete outgoing energy, with its
   *        associated discrete cosines.
   */
  struct DiscreteEnergy {
    double energy;               /**< Discrete outgoing energy */
    std::vector<double> cosines; /**< Discrete cosines */
  };

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  std::optional<double> angle_pdf(double E_in, double mu) const override final;

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override final;

  /**
   * @brief Returns true if the outgoing energies are skewed.
   */
  bool skewed() const { return skewed_; }

  /**
   * @brief Returns vector to the incoming energy grid.
   */
  const std::vector<double>& incoming_energy() const {
    return incoming_energy_;
  }

  /**
   * @brief Returns the vector of outgoing energies for all incoming energies.
   */
  const std::vector<std::vector<DiscreteEnergy>>& outgoing_energies() const {
    return outgoing_energies_;
  }

 private:
  std::vector<double> incoming_energy_;
  std::vector<std::vector<DiscreteEnergy>> outgoing_energies_;
  uint32_t Noe, Nmu;
  bool skewed_;
};

}  // namespace pndl

#endif
