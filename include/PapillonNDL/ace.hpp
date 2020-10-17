/*
 * Copyright 2020, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * Ce  logiciel  est  un  programme  informatique  servant  à  résoudre
 * l'équation de Boltzmann  pour les  neutrons, avec  des énergies continus,
 * en trois  dimensions  spatiales, utilisant la méthode  Monte Carlo.
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

#include <array>
#include <string>
#include <vector>

namespace pndl {
class ACE {
 public:
  ACE(std::string fname);
  ~ACE() = default;

  // Basic nuclide data
  int32_t zaid() const;
  double temperature() const;
  double awr() const;
  bool fissile() const;

  // Accessors to data arrays
  std::pair<int32_t, double> izaw(size_t i) const;
  int32_t nxs(size_t i) const;
  int32_t jxs(size_t i) const;
  double xss(size_t i) const;

  template <class T>
  T xss(size_t i) const {
    return static_cast<T>(xss(i));
  }

  // Range accessors to data arrays
  std::vector<std::pair<int32_t, double>> izaw(size_t i, size_t len) const;
  std::vector<int32_t> nxs(size_t i, size_t len) const;
  std::vector<int32_t> jxs(size_t i, size_t len) const;
  std::vector<double> xss(size_t i, size_t len) const;

  template <class T>
  std::vector<T> xss(size_t i, size_t len) const {
    std::vector<T> tmp(len);

    for (size_t j = 0; j < len; j++) {
      tmp[j] = static_cast<T>(xss<T>(i + j));
    }

    return tmp;
  }

  // Locations for sections within XSS
  int32_t ESZ() const;
  int32_t NU() const;
  int32_t MTR() const;
  int32_t LQR() const;
  int32_t TYR() const;
  int32_t LSIG() const;
  int32_t SIG() const;
  int32_t LAND() const;
  int32_t AND() const;
  int32_t LDLW() const;
  int32_t DLW() const;
  //=======================================================================

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

};  // ACE
}  // namespace pndl

#endif  // PAPILLON_NDL_ACE_H
