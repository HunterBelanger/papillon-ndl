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
#ifndef PAPILLON_NDL_LEVEL_INELASTIC_SCATTER_H
#define PAPILLON_NDL_LEVEL_INELASTIC_SCATTER_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <memory>

namespace pndl {

/**
 * @brief Energy distribution for inelastic scatter.
 */
class LevelInelasticScatter : public EnergyLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   */
  LevelInelasticScatter(const ACE& ace, std::size_t i);

  /**
   * @param Q Q-value of the reaction.
   * @param AWR Atomic Weight Ratio of nuclide.
   */
  LevelInelasticScatter(double Q, double AWR);
  ~LevelInelasticScatter() = default;

  double sample_energy(double E_in,
                       std::function<double()> rng) const override final;

  double pdf(double E_in, double E_out) const override final;

  /**
   * @brief Returns first parameter which is (A+1)*abs(Q)/A.
   */
  double C1() const {return C1_;}

  /**
   * @brief Returns second parameter which is (A/(A+1))^2.
   */
  double C2() const {return C2_;}

 private:
  double C1_;
  double C2_;
};

}  // namespace pndl

#endif
