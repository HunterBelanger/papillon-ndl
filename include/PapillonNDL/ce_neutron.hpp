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
#ifndef PAPILLON_NDL_CE_NEUTRON_H
#define PAPILLON_NDL_CE_NEUTRON_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ce_neutron_base.hpp>
#include <PapillonNDL/reaction.hpp>

namespace pndl {

template <typename XSType>
class CENeutron {};

template <>
class CENeutron<CrossSection> : public CENeutronBase {
 public:
  /**
   * @param ace ACE file from which to construct the data.
   */
  CENeutron(const ACE& ace);

  /**
   * @param ace ACE file from which to take the new cross sections.
   * @param nuclide CENeutron containing another instance of the desired
   *                nuclide. Secondary distributions and fission data
   *                will be shared between the two data sets.
   */
  CENeutron(const ACE& ace, const CENeutron& nuclide);

  /**
   * @brief Returns the temperature at which the data has been prepared.
   */
  double temperature() const { return temperature_; }

  /**
   * @brief Returns the energy grid for the nuclide.
   */
  const EnergyGrid& energy_grid() const { return *energy_grid_; }

  /**
   * @brief Returns the total CrossSection for the nuclide.
   */
  const CrossSection& total_xs() const { return *total_xs_; }

  /**
   * @brief Returns the elastic scattering CrossSection for the
   * nuclide.
   */
  const CrossSection& elastic_xs() const { return *elastic_xs_; }

  /**
   * @brief Returns the fission CrossSection for the nuclide.
   */
  const CrossSection& fission_xs() const { return *fission_xs_; }

  /**
   * @brief Returns the disappearance CrossSection for the nuclide.
   */
  const CrossSection& disappearance_xs() const { return *disappearance_xs_; }

  /**
   * @brief Returns the photon production CrossSection for the nuclide.
   */
  const CrossSection& photon_production_xs() const {
    return *photon_production_xs_;
  }

  /**
   * @brief Retrieved a given MT reaction.
   * @param mt MT reaction to return.
   */
  const STReaction& reaction(uint32_t mt) const {
    if (!this->has_reaction(mt)) {
      std::string mssg = "MT = " + std::to_string(mt) +
                         " is not provided in ZAID = " + std::to_string(zaid_) +
                         ".";
      throw PNDLException(mssg);
    }

    return reactions_[reaction_indices_[mt]];
  }

 private:
  double temperature_;

  std::shared_ptr<EnergyGrid> energy_grid_;
  std::shared_ptr<CrossSection> total_xs_;
  std::shared_ptr<CrossSection> disappearance_xs_;
  std::shared_ptr<CrossSection> elastic_xs_;
  std::shared_ptr<CrossSection> fission_xs_;
  std::shared_ptr<CrossSection> photon_production_xs_;

  std::vector<STReaction> reactions_;

  // Private Helper Methods
  std::shared_ptr<CrossSection> compute_fission_xs();
};

using STNeutron = CENeutron<CrossSection>;

}  // namespace pndl

#endif
