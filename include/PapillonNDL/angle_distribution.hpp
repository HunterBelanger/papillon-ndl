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
#ifndef PAPILLON_NDL_ANGLE_DISTRIBUTION_H
#define PAPILLON_NDL_ANGLE_DISTRIBUTION_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_law.hpp>
#include <functional>
#include <memory>

namespace pndl {

/**
 *  @brief Holds all of the angular distributions at all provided energies
 *         for a single reaction.
 */
class AngleDistribution {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param locb Index in the XSS array to start reading angular distribution.
   */
  AngleDistribution(const ACE& ace, int locb);

  /**
   * @param energy_grid Incoming energies which has angle laws.
   * @param laws Angle laws for each incoming energy.
   */
  AngleDistribution(const std::vector<double>& energy_grid,
                    const std::vector<std::shared_ptr<AngleLaw>>& laws);
  ~AngleDistribution() = default;

  /**
   * @brief Samples a scattering cosine for the given energy.
   * @param E_in Incident energy before scatter, in MeV.
   * @param rng Random number generator function.
   */
  double sample_angle(double E_in, std::function<double()> rng) const;

  /**
   * @brief Evaluates the PDF for having a scattering cosine of mu at incoming
   *        energy E_in.
   * @param E_in Incoming energy.
   * @param mu Scattering cosine.
   */
  double pdf(double E_in, double mu) const;

  /**
   * @brief Returns the number of energies/angular distributions stored.
   */
  std::size_t size() const {return energy_grid_.size();}

  /**
   * @brief Reference to the vector of energy values (in MeV)
   *        which have an angular distribution.
   */
  const std::vector<double>& energy() const {return energy_grid_;}

  /**
   * @brief Gets the ith energy point (in MeV) which has an angular
   * distribution.
   * @param i Index in the energy grid.
   */
  double energy(std::size_t i) const {return energy_grid_[i];}

  /**
   * @brief Gets a pointer to the angular distribution for the ith energy point.
   * @param i Index in the energy grid.
   */
  const AngleLaw& law(std::size_t i) const {return *laws_[i];}

 private:
  std::vector<double> energy_grid_;
  std::vector<std::shared_ptr<AngleLaw>> laws_;
};

}  // namespace pndl

#endif
