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
#ifndef PAPILLON_NDL_ST_INCOHERENT_INELASTIC_H
#define PAPILLON_NDL_ST_INCOHERENT_INELASTIC_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/region_1d.hpp>

namespace pndl {

/**
 * @brief Holds the Incoherent Inelastic scattering data for a single nuclide
 *        at a single temperature.
 */
class STIncoherentInelastic {
 public:
  /**
   * @param ace ACE file which contains thermal scattering law.
   * @param unit_based_interpolation If false (default value), the distribution
   *        will be sampled without using unit-based interpolation, which is
   *        the method used by MCNP, Serpent, and OpenMC. If set to true, unit
   *        based interpolation will be applied to the sampling of the energy.
   */
  STIncoherentInelastic(const ACE& ace, bool unit_based_interpolation = false);
  ~STIncoherentInelastic() = default;

  /**
   * @brief Returns a pointer to the cross section function.
   */
  std::shared_ptr<Region1D> cross_section() const { return xs_; }

  /**
   * @brief Evaluates the incoherent inelastic scattering cross section
   *        at energy E.
   * @param E Incident energy at which to evaluate the cross section in MeV.
   */
  double xs(double E) const { return (*this->xs_)(E); }

  /**
   * @brief Retruns the maximum energy value which is tabulated for the
   *        cross section.
   */
  double max_energy() const { return this->xs_->max_x(); }

  /**
   * @brief Sample the angle-energy distribution.
   * @param E_in Incident energy in MeV.
   * @param rng Random number generation function.
   */
  AngleEnergyPacket sample_angle_energy(double E_in,
                                        std::function<double()> rng) const {
    return angle_energy_->sample_angle_energy(E_in, rng);
  }

  /**
   * @brief Returns a pointer to the AngleEnergy distribution.
   */
  std::shared_ptr<AngleEnergy> distribution() const { return angle_energy_; }

 private:
  std::shared_ptr<Region1D> xs_;
  std::shared_ptr<AngleEnergy> angle_energy_;
};

}  // namespace pndl

#endif
