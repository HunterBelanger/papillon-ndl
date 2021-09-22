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
#ifndef PAPILLON_NDL_CE_NEUTRON_BASE_H
#define PAPILLON_NDL_CE_NEUTRON_BASE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_distribution.hpp>
#include <PapillonNDL/delayed_group.hpp>
#include <array>
#include <memory>

namespace pndl {

/**
 * @brief Holds all non temperature-dependent, continuous energy data for a
 *        single nuclide.
 */
class CENeutronBase {
 public:
  CENeutronBase(const CENeutronBase&) = default;

  /**
   * @brief Returns the nuclide ZAID.
   */
  uint32_t zaid() const { return zaid_; }

  /**
   * @brief Returns the nuclide Atomic Weight Ratio.
   */
  double awr() const { return awr_; }

  /**
   * @brief Returns true if the nuclide is fissile, and false otherwise.
   */
  bool fissile() const { return fissile_; }

  /**
   * @brief Returns the function for total nu.
   */
  const Function1D& nu_total() const { return *nu_total_; }

  /**
   * @brief Returns the function for prompt nu.
   */
  const Function1D& nu_prompt() const { return *nu_prompt_; }

  /**
   * @brief Returns the function for delayed nu.
   */
  const Function1D& nu_delayed() const { return *nu_delayed_; }

  /**
   * @brief Returns the AngleDistribution for elastic scattering.
   */
  const AngleDistribution& elastic_angle_distribution() const {
    return *elastic_angle_;
  }

  /**
   * @brief Returns the number of delayed neutron groups.
   */
  std::size_t n_delayed_groups() const { return delayed_groups_.size(); }

  /**
   * @brief Returns the ith delayed group data.
   * @param i Index of the delayed group.
   */
  const DelayedGroup& delayed_group(std::size_t i) const {
    return delayed_groups_[i];
  }

  /**
   * @brief Returns a list of all MT reactions present for the nuclide.
   */
  const std::vector<uint32_t>& mt_list() const { return mt_list_; }

  /**
   * @brief Checks to see if a nucldie has a given reaction.
   * @param mt MT reaction to search for.
   */
  bool has_reaction(uint32_t mt) const {
    return (mt > 891 || reaction_indices_[mt] < 0) ? false : true;
  }

 protected:
  uint32_t zaid_;
  double awr_;
  bool fissile_;

  std::shared_ptr<AngleDistribution> elastic_angle_;

  std::shared_ptr<Function1D> nu_total_;
  std::shared_ptr<Function1D> nu_prompt_;
  std::shared_ptr<Function1D> nu_delayed_;
  std::vector<DelayedGroup> delayed_groups_;

  std::vector<uint32_t> mt_list_;
  std::array<int32_t, 892> reaction_indices_;

  /**
   * @param ace ACE file from which to construct the data.
   */
  CENeutronBase(const ACE& ace);

  // Private helper methods
  void read_fission_data(const ACE& ace);
  std::shared_ptr<Function1D> read_nu(const ACE& ace, std::size_t i);
  std::shared_ptr<Function1D> read_polynomial_nu(const ACE& ace, std::size_t i);
  std::shared_ptr<Function1D> read_tabular_nu(const ACE& ace, std::size_t i);
};

}  // namespace pndl

#endif
