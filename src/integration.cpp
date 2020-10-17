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
#include <PapillonNDL/integration.hpp>
#include <cmath>
#include <stdexcept>

namespace pndl {

double histogram_integrate(double x_low, double x_hi, double /*x1*/, double y1,
                           double /*x2*/, double /*y2*/) {
  return y1 * (x_hi - x_low);
}

double linear_linear_integrate(double x_low, double x_hi, double x1, double y1,
                               double x2, double y2) {
  // Here, y is linear in x
  // y = ((x - x1)/(x2 - x1))*(y2 - y1) + y1
  double numerator = (x_hi - x_low) * (y1 - y2) * (x_hi + x_low - 2. * x1);
  double denominator = 2. * (x1 - x2);
  return (numerator / denominator) + (x_hi - x_low) * y1;
}

double log_linear_integrate(double x_low, double x_hi, double x1, double y1,
                            double x2, double y2) {
  // Here, y is linear in log(x)
  // log(y) = ((x - x1)/(x2 - x1))*log(y2/y1) + log(y1)
  if (y2 / y1 <= 0.) {
    throw std::runtime_error(
        "Integration: log_linear: Must satisfy y2 / y1 > 0.");
  }

  double base = y2 / y1;
  double denominator = std::log(base);
  double coefficient = y1 * (x2 - x1);
  double exponent_hi = (x1 - x_hi) / (x1 - x2);
  double exponent_low = (x1 - x_low) / (x1 - x2);
  return (coefficient / denominator) *
         (std::pow(base, exponent_hi) - std::pow(base, exponent_low));
}

double linear_log_integrate(double x_low, double x_hi, double x1, double y1,
                            double x2, double y2) {
  // Here, y is linear in log(x)
  // y = (log(x/x1)/log(x2/x1))*(y2 - y1) + y1
  if (x_hi / x1 <= 0.) {
    throw std::runtime_error(
        "Integration: linear_log: Must satisfy x_hi / x1 > 0.");
  }

  if (x_low / x1 <= 0.) {
    throw std::runtime_error(
        "Integration: linear_log: Must satisfy x_low / x1 > 0.");
  }

  if (x2 / x1 <= 0.) {
    throw std::runtime_error(
        "Integration: linear_log: Must satisfy x2 / x1 > 0.");
  }

  double numerator_hi = x_hi * ((y2 - y1) * std::log(x_hi / x1) +
                                y1 * std::log(x2 / x1) + y1 - y2);
  double numerator_low = x_low * ((y2 - y1) * std::log(x_low / x1) +
                                  y1 * std::log(x2 / x1) + y1 - y2);
  double denominator = std::log(x2 / x1);
  return (numerator_hi / denominator) - (numerator_low / denominator);
}

double log_log_integrate(double x_low, double x_hi, double x1, double y1,
                         double x2, double y2) {
  if (y2 / y1 <= 0.) {
    throw std::runtime_error("Integration: log_log: Must satisfy y2 / y1 > 0.");
  }

  if (x2 / x1 <= 0.) {
    throw std::runtime_error("Integration: log_log: Must satisfy x2 / x1 > 0.");
  }

  double y2_y1 = y2 / y1;
  double x2_x1 = x2 / x1;
  double log_y2_y1 = std::log(y2_y1);
  double log_x2_x1 = std::log(x2_x1);
  double exponent = log_y2_y1 / log_x2_x1;
  double denominator = exponent + 1.0;
  return (y1 / denominator) * (x_hi * std::pow(x_hi / x1, exponent) -
                               x_low * std::pow(x_low / x1, exponent));
}

double integrate(double x_low, double x_hi, double x1, double y1, double x2,
                 double y2, Interpolation interp) {
  if (x_low < x1 || x_low > x2) {
    throw std::runtime_error("Integration: Must satisfy x1 <= x_low <= x2");
  }

  if (x_hi < x1 || x_hi > x2) {
    throw std::runtime_error("Integration: Must satisfy x1 <= x_hi <= x2");
  }

  switch (interp) {
    case Interpolation::Histogram:
      return histogram_integrate(x_low, x_hi, x1, y1, x2, y2);
    case Interpolation::LinLin:
      return linear_linear_integrate(x_low, x_hi, x1, y1, x2, y2);
    case Interpolation::LinLog:
      return linear_log_integrate(x_low, x_hi, x1, y1, x2, y2);
    case Interpolation::LogLin:
      return log_linear_integrate(x_low, x_hi, x1, y1, x2, y2);
    case Interpolation::LogLog:
      return log_log_integrate(x_low, x_hi, x1, y1, x2, y2);
  }

  // NEVER REACHED
  return 0;
}

}  // namespace pndl
