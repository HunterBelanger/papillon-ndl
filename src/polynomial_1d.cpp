/*
 * Copyright 2020, Hunter Belanger
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
#include <PapillonNDL/polynomial_1d.hpp>

namespace pndl {

Polynomial1D::Polynomial1D(const std::vector<double>& coeffs)
    : coefficients_(coeffs) {}

double Polynomial1D::operator()(double x) const {
  double value = 0.;
  for (auto coeff = coefficients_.rbegin(); coeff != coefficients_.rend();
       coeff++) {
    value = (value * x) + *coeff;
  }
  return value;
}

double Polynomial1D::integrate(double x_low, double x_hi) const {
  double integral_hi = 0.;
  double integral_low = 0.;

  double exponent = static_cast<double>(order() + 1);
  for (auto coeff_it = coefficients_.rbegin(); coeff_it != coefficients_.rend();
       coeff_it++) {
    double coeff = *coeff_it / exponent;
    integral_hi = (integral_hi * x_hi) + coeff;
    integral_low = (integral_low * x_low) + coeff;
    exponent -= 1.;
  }
  // Do one extra where coeff = 0
  integral_hi = integral_hi * x_hi;
  integral_low = integral_low * x_low;

  return integral_hi - integral_low;
}

size_t Polynomial1D::order() const { return coefficients_.size() - 1; }

double Polynomial1D::coefficient(size_t i) const { return coefficients_[i]; }

}  // namespace pndl
