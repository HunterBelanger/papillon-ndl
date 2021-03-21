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
#ifndef PAPILLON_NDL_ANGLE_TABLE_H
#define PAPILLON_NDL_ANGLE_TABLE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/angle_law.hpp>
#include <PapillonNDL/pctable.hpp>

namespace pndl {

/**
 * @brief Angular distribution which is provided as tabulated PDF and CDF.
 */
class AngleTable : public AngleLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   */
  AngleTable(const ACE& ace, size_t i);

  /**
   * @param cosines Conines of scattering angle which are tabulated.
   * @param pdf The Probability Density Function for the provided values.
   * @param cdf The Cumulative Density Function for the provided values.
   * @param interp Interpolation rule for the data. May be either
   *               Histogram or LinLin.
   */
  AngleTable(const std::vector<double>& cosines, const std::vector<double>& pdf,
             const std::vector<double>& cdf, Interpolation interp);

  /**
   * @param table PCTable contianing the PDF and CDF for the cosine
   *              distribution.
   */
  AngleTable(const PCTable& table);
  ~AngleTable() = default;

  double sample_mu(double xi) const override final;

  /**
   * @brief Returns the number of points in the tabulated data.
   */
  size_t size() const;

  /**
   * @brief Returns the vector of the cosine points.
   */
  const std::vector<double>& cosines() const;

  /**
   * @brief Returns the vector of the PDF values.
   */
  const std::vector<double>& pdf() const;

  /**
   * @brief Returns the vector of the CDF values.
   */
  const std::vector<double>& cdf() const;

  /**
   * @brief Returns the type of interpolation used on the table
   *        (Histogram or LinLin).
   */
  Interpolation interpolation() const;

 private:
  PCTable distribution_;
};

}  // namespace pndl

#endif
