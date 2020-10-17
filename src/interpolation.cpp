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
#include <PapillonNDL/interpolation.hpp>
#include <cmath>
#include <stdexcept>

namespace pndl {

double histogram(double x, double x1, double y1, double x2, double y2) {
  if (x >= x1 && x < x2)
    return y1;
  else
    return y2;
}

double linear_linear(double x, double x1, double y1, double x2, double y2) {
  // Here, y is linear in x
  // y = ((x - x1)/(x2 - x1))*(y2 - y1) + y1
  double f = interpolation_factor(x, x1, x2);
  return (y2 - y1) * f + y1;
}

double log_linear(double x, double x1, double y1, double x2, double y2) {
  // Here, y is linear in log(x)
  // log(y) = ((x - x1)/(x2 - x1))*log(y2/y1) + log(y1)
  return std::exp(linear_linear(x, x1, std::log(y1), x2, std::log(y2)));
}

double linear_log(double x, double x1, double y1, double x2, double y2) {
  // Here, y is linear in log(x)
  // y = (log(x/x1)/log(x2/x1))*(y2 - y1) + y1
  return linear_linear(std::log(x), std::log(x1), y1, std::log(x2), y2);
}

double log_log(double x, double x1, double y1, double x2, double y2) {
  return std::exp(linear_linear(std::log(x), std::log(x1), std::log(y1),
                                std::log(x2), std::log(y2)));
}

double interpolation_factor(double x, double x1, double x2) {
  return (x - x1) / (x2 - x1);
}

double interpolate(double x, double x1, double y1, double x2, double y2,
                   Interpolation interp) {
  if (x < x1) {
    throw std::runtime_error("Interpolation: Must satisfy x >= x1");
  } else if (x > x2) {
    throw std::runtime_error("Interpolation: Must satisfy x <= x2");
  }

  switch (interp) {
    case Interpolation::Histogram:
      return histogram(x, x1, y1, x2, y2);
    case Interpolation::LinLin:
      return linear_linear(x, x1, y1, x2, y2);
    case Interpolation::LinLog:
      return linear_log(x, x1, y1, x2, y2);
    case Interpolation::LogLin:
      return log_linear(x, x1, y1, x2, y2);
    case Interpolation::LogLog:
      return log_log(x, x1, y1, x2, y2);
  }

  // NEVER REACHED
  return 0;
}

}  // namespace pndl
