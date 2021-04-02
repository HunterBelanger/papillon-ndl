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
#ifndef PAPILLON_NDL_NBODY_H
#define PAPILLON_NDL_NBODY_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>

namespace pndl {

/**
 * @brief Implements product Angle-Energy law which follows an N-body phase
 *        space distribution.
 */
class NBody : public AngleEnergy {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   * @param iQ Q-value for the reaction.
   * @param probability Function for probability of validity.
   */
  NBody(const ACE& ace, size_t i, double iQ,
        std::shared_ptr<Tabulated1D> probability);

  /**
   * @param n Number of particles (3, 4, or 5).
   * @param Ap Total mass ratio for the n particles.
   * @param AWR Atomic Weight Ratio of the nuclide.
   * @param Q Q-value for the reaction.
   * @param probability Function for probability of validity.
   */
  NBody(uint16_t n, double Ap, double AWR, double Q,
        std::shared_ptr<Tabulated1D> probability);
  ~NBody() = default;

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  /**
   * @brief Returns the number of bodies.
   */
  uint32_t n() const;

  /**
   * @brief Returns the total AWR for all of the particles.
   */
  double Ap() const;

  /**
   * @brief Returns the AWR of the nuclide in question.
   */
  double A() const;

  /**
   * @brief Returns the Q-value for the reaction in MeV.
   */
  double Q() const;

 private:
  uint32_t n_;
  double Ap_;
  double A_;
  double Q_;

  double maxwellian_spectrum(std::function<double()>& rng) const;
};

}  // namespace pndl

#endif
