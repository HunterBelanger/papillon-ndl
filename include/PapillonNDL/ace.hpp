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
#ifndef PAPILLON_NDL_ACE_H
#define PAPILLON_NDL_ACE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <array>
#include <string>
#include <vector>

namespace pndl {

/**
 * @brief Contains data from A Compact ENDF file.
 *
 */
class ACE {
 public:
  /**
   * @param fname Name of the file to be loaded.
   */
  ACE(std::string fname);
  ~ACE() = default;

  /**
   * @brief Gets the ZAID of nuclide represented.
   */
  int32_t zaid() const;

  /**
   * @brief Gets the temperature for which the data was prepared, in kelvins.
   */
  double temperature() const;

  /**
   *  @brief Gets the Atomic Weight Ratio (AWR) of the nuclide.
   */
  double awr() const;

  /**
   * @brief Returns true for a fissile nuclide, false for a non-fissile nuclide.
   */
  bool fissile() const;

  /**
   * @brief Retrieves an (int32_t, double) pair from the IZAW array.
   * @param i index to a pair in the IZAW array.
   *          Must be in the range [0,16).
   */
  std::pair<int32_t, double> izaw(size_t i) const;

  /**
   * @brief Retrieves a value from the NXS array.
   * @param i index to element in the NXS array.
   *          Must be in the range [0,16).
   */
  int32_t nxs(size_t i) const;

  /**
   * @brief Retrieves a value from the JXS array.
   * @param i index to element in the JXS array.
   *          Must be in the range [0,32).
   */
  int32_t jxs(size_t i) const;

  /**
   * @brief Retrieves a value from the XSS array as a double.
   * @param i index to element in the XSS array.
   */
  double xss(size_t i) const;

  /**
   * @brief Retrieves a value from the XSS array, cast to type T.
   * @param i index to element in the XSS array.
   */
  template <class T>
  T xss(size_t i) const {
    return static_cast<T>(xss_[i]);
  }

  /**
   * @brief Retrieves a vector contianing a continuous segment of
   *        (int32_t,double) pairs from the IZAW array.
   * @param i Starting index in the IZAW array.
   * @param len Number of elements to return.
   */
  std::vector<std::pair<int32_t, double>> izaw(size_t i, size_t len) const;

  /**
   * @brief Retrieves a vector contianing a continuous segment of
   *        values from the NXS array.
   * @param i Starting index in the NXS array.
   * @param len Number of elements to return.
   */
  std::vector<int32_t> nxs(size_t i, size_t len) const;

  /**
   * @brief Retrieves a vector contianing a continuous segment of
   *        values from the JXS array.
   * @param i Starting index in the JXS array.
   * @param len Number of elements to return.
   */
  std::vector<int32_t> jxs(size_t i, size_t len) const;

  /**
   * @brief Retrieves a vector contianing a continuous segment of
   *        values from the XSS array.
   * @param i Starting index in the XSS array.
   * @param len Number of elements to return.
   */
  std::vector<double> xss(size_t i, size_t len) const;

  /**
   * @brief Retrieves a vector contianing a continuous segment of
   *        values from the XSS array, all cast to type T.
   * @param i Starting index in the XSS array.
   * @param len Number of elements to return.
   */
  template <class T>
  std::vector<T> xss(size_t i, size_t len) const {
    std::vector<T> tmp(len);

    for (size_t j = 0; j < len; j++) {
      tmp[j] = static_cast<T>(xss<T>(i + j));
    }

    return tmp;
  }

  /**
   * @brief Returns a pointer to the beginning of the XSS array.
   */
  const double* xss_data() const;

  /**
   * @brief Returns the index to the beginning of the ESZ block.
   */
  int32_t ESZ() const;

  /**
   * @brief Returns the index to the beginning of the NU block.
   */
  int32_t NU() const;

  /**
   * @brief Returns the index to the beginning of the MTR block.
   */
  int32_t MTR() const;

  /**
   * @brief Returns the index to the beginning of the LQR block.
   */
  int32_t LQR() const;

  /**
   * @brief Returns the index to the beginning of the TYR block.
   */
  int32_t TYR() const;

  /**
   * @brief Returns the index to the beginning of the LSIG block.
   */
  int32_t LSIG() const;

  /**
   * @brief Returns the index to the beginning of the SIG block.
   */
  int32_t SIG() const;

  /**
   * @brief Returns the index to the beginning of the LAND block.
   */
  int32_t LAND() const;

  /**
   * @brief Returns the index to the beginning of the AND block.
   */
  int32_t AND() const;

  /**
   * @brief Returns the index to the beginning of the LDLW block.
   */
  int32_t LDLW() const;

  /**
   * @brief Returns the index to the beginning of the DLW block.
   */
  int32_t DLW() const;

  /**
   * @brief Returns the index to the beginning of the DNEDL block.
   */
  int32_t DNEDL() const;

  /**
   * @brief Returns the index to the beginning of the DNED block.
   */
  int32_t DNED() const;

  /**
   * @brief Returns the index to the beginning of the DNU block.
   */
  int32_t DNU() const;

  /**
   * @brief Returns the index to the beginning of the BDD block.
   */
  int32_t BDD() const;

  /**
   * @brief Returns the index to the beginning of the GPD block.
   */
  int32_t GPD() const;

 private:
  int32_t zaid_;
  double temperature_;
  double awr_;
  bool fissile_;

  std::array<std::pair<int32_t, double>, 16> izaw_;
  std::array<int32_t, 16> nxs_;
  std::array<int32_t, 32> jxs_;
  std::vector<double> xss_;

  // Locator constants
  int32_t esz_, nu_, mtr_, lqr_, tyr_, lsig_, sig_, lan_, an_, ldlw_, dlw_;
  int32_t dnedl_, dned_, dnu_, bdd_, gpd_;

};  // ACE
}  // namespace pndl

#endif  // PAPILLON_NDL_ACE_H
