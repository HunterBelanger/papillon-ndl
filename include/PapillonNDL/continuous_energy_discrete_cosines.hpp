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
   */
  ContinuousEnergyDiscreteCosines(const ACE& ace);

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
  };

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  /**
   * @brief Returns vector to the incoming energy grid.
   */
  const std::vector<double>& incoming_energy() const {
    return incoming_energy_;
  }

  /**
   * @brief Returns the number of incoming energy points.
   */
  size_t size() const { return incoming_energy_.size(); }

  /**
   * @brief Returns the vector of all CEDCTables.
   */
  const std::vector<CEDCTable> tables() const { return tables_; }

  /**
   * @brief Returns a CEDCTable which contains the distributions
   *        for the ith incoming energy.
   * @param i Index to the incoming energy.
   */
  const CEDCTable& table(size_t i) const { return tables_[i]; }

 private:
  std::vector<double> incoming_energy_;
  std::vector<CEDCTable> tables_;
  uint32_t Nmu;

  // Table to use for determining outgoing energy, and j
  // is the index into the cosine array, for sampling the
  // cosine of the scattering angle.
  double sample_energy(const CEDCTable& table, double xi, size_t& j) const {
    double E_out = 0.;
    auto cdf_it = std::lower_bound(table.cdf.begin(), table.cdf.end(), xi);
    if (cdf_it == table.cdf.begin()) {
      E_out = table.energy.front();
      j = 0;
      return E_out;
    }
    size_t l = std::distance(table.cdf.begin(), cdf_it) - 1;

    l = std::min(l, table.energy.size() - 2);

    // Must account for case where pdf_[l] = pdf_[l+1], which means  that
    // the slope is zero, and m=0. This results in nan for the linear alg.
    // To avoid this, must use histogram for that segment.
    if (table.pdf[l] == table.pdf[l + 1]) {
      E_out = histogram_interp_energy(table, xi, l);
    } else {
      E_out = linear_interp_energy(table, xi, l);
    }

    // Set j to be l, so we know which cosines to use latter
    // when sampling the angle.
    j = l;

    return E_out;
  }

  double histogram_interp_energy(const CEDCTable& table, double xi,
                                 size_t l) const {
    return table.energy[l] + ((xi - table.cdf[l]) / table.pdf[l]);
  }

  double linear_interp_energy(const CEDCTable& table, double xi,
                              size_t l) const {
    double m = (table.pdf[l + 1] - table.pdf[l]) /
               (table.energy[l + 1] - table.energy[l]);

    return table.energy[l] +
           (1. / m) * (std::sqrt(table.pdf[l] * table.pdf[l] +
                                 2. * m * (xi - table.cdf[l])) -
                       table.pdf[l]);
  }
};

}  // namespace pndl

#endif