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
#ifndef PAPILLON_NDL_EQUIPROBABLE_ENERGY_BINS_H
#define PAPILLON_NDL_EQUIPROBABLE_ENERGY_BINS_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_law.hpp>

namespace pndl {

/**
 * @brief Energy distribution which is provided as equiprobable energy bins.
 */
class EquiprobableEnergyBins : public EnergyLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   */
  EquiprobableEnergyBins(const ACE& ace, size_t i);
  ~EquiprobableEnergyBins() = default;

  double sample_energy(double E_in,
                       std::function<double()> rng) const override final;

  /**
   * @brief Returns a vector of the grid of incoming energy points for which
   *        an equiprobable bin set is stored.
   */
  const std::vector<double>& incoming_energy() const;

  /**
   * @brief Returns the ith set of bin boundaries as a vector.
   * @param i Index for the incoming energy grid.
   */
  const std::vector<double>& bin_bounds(size_t i) const;

  /**
   * @brief Returns the number of incoming energies / bin boundary sets stored.
   */
  size_t size() const;

 private:
  std::vector<double> incoming_energy_;
  std::vector<std::vector<double>> bin_sets_;

  double sample_bins(double xi1, double xi2,
                     const std::vector<double>& bounds) const;
};

}  // namespace pndl

#endif
