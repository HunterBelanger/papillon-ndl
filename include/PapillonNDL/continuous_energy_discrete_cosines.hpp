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
#ifndef PAPILLON_NDL_CONTINUOUS_ENERGY_DISCRETE_COSINES_H
#define PAPILLON_NDL_CONTINUOUS_ENERGY_DISCRETE_COSINES_H

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>

namespace pndl {

/**
 * @brief Class which represents continuous energy distributions with
 *        discrete cosines for incoherent inelastic scattering.
 */
class ContinuousEnergyDiscreteCosines : public AngleEnergy {
 public:
  /**
   * @param ace ACE file which contains thermal scattering law.
   * @param unit_based_interpolation If false (default value), the distribution
   *        will be sampled without using unit-based interpolation, which is
   *        the method used by MCNP, Serpent, and OpenMC. If set to true, unit
   *        based interpolation will be applied to the sampling of the energy.
   */
  ContinuousEnergyDiscreteCosines(const ACE& ace,
                                  bool unit_based_interpolation = false);

  /**
   * @brief Struct which contains the outgoing energy distribution and
   *        the discrete scattering cosiens for a single incident energy.
   */
  struct CEDCTable {
    std::vector<double> energy; /**< Outgoing energy points */
    std::vector<double> pdf;    /**< PDF for the outgoing energy */
    std::vector<double> cdf;    /**< CDF for the outgoing energy */
    std::vector<std::vector<double>>
        cosines; /**< Discrete scattering cosines for each outgoing energy */

    /**
     * @brief Samples and outgoing energy from the distributions while
     *        also setting the value j to be the lower bound index
     *        for sampling the scattering cosine latter on.
     * @param xi Random variable on the unit interval [0,1).
     * @param j Index which will be set to latter locate the proper
     *          distribution for the scattering cosine.
     */
    double sample_energy(double xi, std::size_t& j) const;
  };

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  std::optional<double> angle_pdf(double E_in, double mu) const override final;

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override final;

  /**
   * @brief Returns vector to the incoming energy grid.
   */
  const std::vector<double>& incoming_energy() const {
    return incoming_energy_;
  }

  /**
   * @brief Returns the number of incoming energy points.
   */
  std::size_t size() const { return incoming_energy_.size(); }

  /**
   * @brief Returns the vector of all CEDCTables.
   */
  const std::vector<CEDCTable> tables() const { return tables_; }

  /**
   * @brief Returns a CEDCTable which contains the distributions
   *        for the ith incoming energy.
   * @param i Index to the incoming energy.
   */
  const CEDCTable& table(std::size_t i) const { return tables_[i]; }

  /**
   * @brief Returns true if the distribution used unit-based interpolation
   *        in sampling the scattering energy and angle, and false otherwise.
   */
  bool unit_based_interpolation() const { return unit_based_interpolation_; }

 private:
  std::vector<double> incoming_energy_;
  std::vector<CEDCTable> tables_;
  uint32_t Nmu;
  bool unit_based_interpolation_;

  AngleEnergyPacket sample_with_unit_based_interpolation(
      double E_in, std::function<double()> rng) const;
  AngleEnergyPacket sample_without_unit_based_interpolation(
      double E_in, std::function<double()> rng) const;
};

}  // namespace pndl

#endif
