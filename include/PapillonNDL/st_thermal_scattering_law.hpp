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
#ifndef PAPILLON_NDL_ST_THERMAL_SCATTERING_LAW_H
#define PAPILLON_NDL_ST_THERMAL_SCATTERING_LAW_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/st_coherent_elastic.hpp>
#include <PapillonNDL/st_incoherent_elastic.hpp>
#include <PapillonNDL/st_incoherent_inelastic.hpp>

namespace pndl {

/**
 * @brief Class to hold all thermal scattering data for a single nuclide
 *        at at single temperature.
 */
class STThermalScatteringLaw {
 public:
  /**
   * @param ace ACE file which contains thermal scattering law.
   */
  STThermalScatteringLaw(const ACE& ace);
  ~STThermalScatteringLaw() = default;

  /**
   * @brief Returns the nuclide ZAID.
   */
  uint32_t zaid() const { return zaid_; }

  /**
   * @brief Returns the nuclide Atomic Weight Ratio.
   */
  double awr() const { return awr_; }

  /**
   * @brief Returns the temperature at which the data has been prepared.
   */
  double temperature() const { return temperature_; }

  /**
   * @brief Returns true if the nuclide has coherrent elastic scattering.
   */
  bool has_coherent_elastic() const {
    if (coherent_elastic_) return true;
    return false;
  }

  /**
   * @brief Returns true if the nuclide has incoherrent elastic scattering.
   */
  bool has_incoherent_elastic() const {
    if (incoherent_elastic_) return true;
    return false;
  }

  /**
   * @brief Returns pointer to the STCoherentElastic instance.
   */
  std::shared_ptr<STCoherentElastic> coherent_elastic() const {
    return coherent_elastic_;
  }

  /**
   * @brief Returns pointer to the STIncoherentElastic instance.
   */
  std::shared_ptr<STIncoherentElastic> incoherent_elastic() const {
    return incoherent_elastic_;
  }

  /**
   * @brief Returns pointer to the STIncoherentInelastic instance.
   */
  std::shared_ptr<STIncoherentInelastic> incoherent_inelastic() const {
    return incoherent_inelastic_;
  }

  /**
   * @brief Returns the value of the coherent elastic scattering
   *        cross section.
   * @param E Energy at which to evaluate the cross section.
   */
  double coherent_elastic_xs(double E) const {
    if (coherent_elastic_) return coherent_elastic_->xs(E);
    return 0.;
  }

  /**
   * @brief Returns the value of the incoherent elastic scattering
   *        cross section.
   * @param E Energy at which to evaluate the cross section.
   */
  double incoherent_elastic_xs(double E) const {
    if (incoherent_elastic_) return incoherent_elastic_->xs(E);
    return 0.;
  }

  /**
   * @brief Returns the value of the incoherent inelastic scattering
   *        cross section.
   * @param E Energy at which to evaluate the cross section.
   */
  double incoherent_inelastic_xs(double E) const {
    if (incoherent_inelastic_) return incoherent_inelastic_->xs(E);
    return 0.;
  }

  /**
   * @breif Returns the total thermal scattering cross section.
   * @param E Energy at which to evaluate the cross section.
   */
  double xs(double E) const {
    return incoherent_inelastic_xs(E) + incoherent_elastic_xs(E) +
           coherent_elastic_xs(E);
  }

 private:
  uint32_t zaid_;
  double awr_;
  double temperature_;

  std::shared_ptr<STCoherentElastic> coherent_elastic_;
  std::shared_ptr<STIncoherentElastic> incoherent_elastic_;
  std::shared_ptr<STIncoherentInelastic> incoherent_inelastic_;
};
}  // namespace pndl

#endif
