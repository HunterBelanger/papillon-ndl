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
#ifndef PAPILLON_NDL_TABULAR_ENERGY_ANGLE_H
#define PAPILLON_NDL_TABULAR_ENERGY_ANGLE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/energy_angle_table.hpp>

namespace pndl {

/**
 * @brief A product Angle-Energy distribution where the angle and energy
 *        PDFs and CDFs are tabulated for different incoming energies.
 */
class TabularEnergyAngle : public AngleEnergy {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   * @param JED Relative index for finding energy and angle distributions.
   */
  TabularEnergyAngle(const ACE& ace, size_t i, size_t JED);

  /**
   * @param incoming_energy Incoming energy grid.
   * @param tables vector of EnergyAngleTable for each point in the
   *               incoming energy grid.
   */
  TabularEnergyAngle(const std::vector<double>& incoming_energy,
                     const std::vector<EnergyAngleTable>& tables);
  ~TabularEnergyAngle() = default;

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  /**
   * @brief Returns a vector to the grid of incoming energies.
   */
  const std::vector<double>& incoming_energy() const;

  /**
   * @brief Returns the ith incoming energy in MeV.
   * @param i Index to the incoming energy grid.
   */
  double incoming_energy(size_t i) const;

  /**
   * @brief Returns an EnergyAngleTable which contains the distributions
   *        for the ith incoming energy.
   * @param i Index to the incoming energy.
   */
  const EnergyAngleTable& table(size_t i) const;

  /**
   * @brief Returns the number of incoming energy points.
   */
  size_t size() const;

 private:
  std::vector<double> incoming_energy_;
  std::vector<EnergyAngleTable> tables_;
};

}  // namespace pndl

#endif
