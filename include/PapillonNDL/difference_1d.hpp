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
#ifndef PAPILLON_NDL_DIFFERENCE_1D_H
#define PAPILLON_NDL_DIFFERENCE_1D_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/function_1d.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <memory>

namespace pndl {

/**
 * @brief Class to represent functions which is the difference of two
 *        other functions.
 */
class Difference1D : public Function1D {
 public:
  /**
   * @brief Function will be evaluated as term1(x) - term2(x).
   * @param term1 Pointer to the function for the first term.
   * @param term2 Pointer to the function for the second term.
   */
  Difference1D(std::shared_ptr<Function1D> term1, std::shared_ptr<Function1D> term2): term_1_(term1), term_2_(term2) {
    if (!term_1_) {
      std::string mssg = "Difference1D::Difference1D: Term 1 is nullptr.";
      throw PNDLException(mssg, __FILE__, __LINE__); 
    }

    if (!term_2_) {
      std::string mssg = "Difference1D::Difference1D: Term 2 is nullptr.";
      throw PNDLException(mssg, __FILE__, __LINE__); 
    }
  }

  double operator()(double x) const override final {
    return (*term_1_)(x) - (*term_2_)(x); 
  }

  double integrate(double x_low, double x_hi) const override final {
    return term_1_->integrate(x_low, x_hi) - term_2_->integrate(x_low, x_hi); 
  }
  
  /**
   * @brief Returns the first function in the difference.
   */
  const Function1D& term_1() const { return *term_1_; }

  /**
   * @brief Returns the second function in the difference.
   */
  const Function1D& term_2() const { return *term_2_; }

 private:
  std::shared_ptr<Function1D> term_1_;
  std::shared_ptr<Function1D> term_2_;
};

}  // namespace pndl

#endif
