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
#ifndef PAPILLON_NDL_INTERPOLATION_H
#define PAPILLON_NDL_INTERPOLATION_H

#include <PapillonNDL/pndl_exception.hpp>
#include <cmath>
#include <cstdint>
#include <algorithm>

namespace pndl {

//==============================================================================
// Interpolation Identifiers
enum class Interpolation : uint32_t {
  Histogram = 1,
  LinLin = 2,
  LinLog = 3,
  LogLin = 4,
  LogLog = 5
};

//==============================================================================
// Helper function to test if values in an interval change sign. Used to
// validate grids for LinLog, LogLin, and LogLog interpolations.
template<class ForwardIt>
bool has_sign_change(ForwardIt first, ForwardIt last) {
  // Get initial sign of array
  auto first_sign = std::signbit(*first);
  for(auto it = first; it != last; it++) {
    if(std::signbit(*it) != first_sign) {
      // There is a sign change in the grid !
      return true;
    }
  }
  return false;
}

//==============================================================================
// Histogram Interpolation
struct Histogram {
  template<class T>
  static T interpolate(T /*x*/, T /*x1*/, T y1, T /*x2*/, T /*y2*/) {
    return y1;
  }

  template<class T>
  static T invert(T /*y*/, T x1, T /*y1*/, T /*x2*/, T /*y2*/) {
    return x1;
  }

  template<class T>
  static T integrate(T x_low, T x_hi, T /*x1*/, T y1, T /*x2*/, T /*y2*/) {
    return y1 * (x_hi - x_low);
  }

  template<class ForwardIt>
  static void verify_x_grid(ForwardIt first, ForwardIt last) {
    auto it = std::adjacent_find(first, last);
    if(it != last) {
      // Grid has repeated elements which are not allowed !
      auto ind = std::distance(first, it);
      std::string mssg = "Histogram::verify_x_grid: Repeated values found in x-grid values\n";
      mssg +=            "of Histogram interpolation at index " + std::to_string(ind) + ".";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  template<class ForwardIt>
  static void verify_y_grid(ForwardIt /*first*/, ForwardIt /*last*/) {
    // No requirements
  }

  static const Interpolation interpolation = Interpolation::Histogram;
};

//==============================================================================
// LinLin Interpolation
struct LinLin {
  template<class T>
  static T interpolate(T x, T x1, T y1, T x2, T y2) {
    return (x - x1)/(x2 - x1) * (y2 - y1) + y1;
  }

  template<class T>
  static T invert(T y, T x1, T y1, T x2, T y2) {
    return (y - y1)/(y2 - y1) * (x2 - x1) + x1;
  }

  template<class T>
  static T integrate(T x_low, T x_hi, T x1, T y1, T x2, T y2) {
    const auto numerator = (x_hi - x_low) * (y1 - y2) * (x_hi + x_low - 2. * x1);
    const auto denominator = 2. * (x1 - x2);
    return (numerator / denominator) + (x_hi - x_low) * y1;
  }

  template<class ForwardIt>
  static void verify_x_grid(ForwardIt first, ForwardIt last) {
    if(!std::is_sorted(first, last)) {
      // Grid elements aren't increasing, which isn't allowed for LinLin
      std::string mssg = "LinLin::verify_x_grid: Non-increasing x-grid values in LinLin\n";
      mssg +=            "interpolation.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  template<class ForwardIt>
  static void verify_y_grid(ForwardIt /*first*/, ForwardIt /*last*/) {
    // No requirements
  }

  static const Interpolation interpolation = Interpolation::LinLin;
};

//==============================================================================
// LinLog Interpolation
struct LinLog {
  template<class T>
  static T interpolate(T x, T x1, T y1, T x2, T y2) {
    return y1 + (y2 - y1) * std::log(x/x1) / std::log(x2/x1);
  }

  template<class T>
  static T invert(T y, T x1, T y1, T x2, T y2) {
    return (y2 != y1) ? x1*std::pow((x2/x1),((y - y1)/(y2 - y1))) : x1;
  }

  template<class T>
  static T integrate(T x_low, T x_hi, T x1, T y1, T x2, T y2) {
    const auto numerator_hi = x_hi * ((y2 - y1) * std::log(x_hi / x1) +
                              y1 * std::log(x2 / x1) + y1 - y2);
    const auto numerator_low = x_low * ((y2 - y1) * std::log(x_low / x1) +
                               y1 * std::log(x2 / x1) + y1 - y2);
    const auto denominator = std::log(x2 / x1);
    return (numerator_hi / denominator) - (numerator_low / denominator);
  }

  template<class ForwardIt>
  static void verify_x_grid(ForwardIt first, ForwardIt last) {
    // Make sure ther are no sign changes
    if(has_sign_change(first, last)) {
      std::string mssg = "LinLog::verify_x_grid: Sign change occurs in x-grid values of LinLog\n";
      mssg +=            "interpolation.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  template<class ForwardIt>
  static void verify_y_grid(ForwardIt /*first*/, ForwardIt /*last*/) {
    // No requirements
  }

  static const Interpolation interpolation = Interpolation::LinLog;
};

//==============================================================================
// LogLin Interpolation
struct LogLin {
  template<class T>
  static T interpolate(T x, T x1, T y1, T x2, T y2) {
    return y1 * std::pow((y2/y1), (x - x1)/(x2 - x1));
  }

  template<class T>
  static T invert(T y, T x1, T y1, T x2, T y2) {
    return (y1 != y2) ? x1+(x2-x1)*std::log(y/y1)/std::log(y2/y1) : x1;
  }

  template<class T>
  static T integrate(T x_low, T x_hi, T x1, T y1, T x2, T y2) {
    const auto base = y2 / y1;
    const auto denominator = std::log(base);
    const auto coefficient = y1 * (x2 - x1);
    const auto exponent_hi = (x1 - x_hi) / (x1 - x2);
    const auto exponent_low = (x1 - x_low) / (x1 - x2);
    return (coefficient / denominator) * 
           (std::pow(base, exponent_hi) - std::pow(base, exponent_low));
  }

  template<class ForwardIt>
  static void verify_x_grid(ForwardIt /*first*/, ForwardIt /*last*/) {
    // No requirements
  }

  template<class ForwardIt>
  static void verify_y_grid(ForwardIt first, ForwardIt last) {
    // Make sure ther are no sign changes
    if(has_sign_change(first, last)) {
      std::string mssg = "LogLin::verify_y_grid: Sign change occurs in y-grid values of LogLin\n";
      mssg +=            "interpolation.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  static const Interpolation interpolation = Interpolation::LogLin;
};

//==============================================================================
// LogLog Interpolation
struct LogLog {
  template<class T>
  static T interpolate(T x, T x1, T y1, T x2, T y2) {
    const auto exponent = std::log(y2/y1) / std::log(x2/x1);
    return y1*std::pow(x/x1, exponent);
  }

  template<class T>
  static T invert(T y, T x1, T y1, T x2, T y2) {
    const auto exponent = std::log(y/y1)/std::log(y2/y1);
    return x1*std::pow(x2/x1, exponent);
  }

  template<class T>
  static T integrate(T x_low, T x_hi, T x1, T y1, T x2, T y2) {
    const auto y2_y1 = y2 / y1;
    const auto x2_x1 = x2 / x1;
    const auto log_y2_y1 = std::log(y2_y1);
    const auto log_x2_x1 = std::log(x2_x1);
    const auto exponent = log_y2_y1 / log_x2_x1;
    const auto denominator = exponent + 1.0;
    if(std::abs(denominator) <= 1.E-12) {return y1*x1*std::log(x_hi/x_low);}
    return (y1 / denominator) * (x_hi * std::pow(x_hi / x1, exponent) -
                                x_low * std::pow(x_low / x1, exponent));
  }

  template<class ForwardIt>
  static void verify_x_grid(ForwardIt first, ForwardIt last) {
    // Make sure ther are no sign changes
    if(has_sign_change(first, last)) {
      std::string mssg = "LogLog::verify_x_grid: Sign change occurs in x-grid values of LogLog\n";
      mssg +=            "interpolation.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  template<class ForwardIt>
  static void verify_y_grid(ForwardIt first, ForwardIt last) {
    // Make sure ther are no sign changes
    if(has_sign_change(first, last)) {
      std::string mssg = "LogLog::verify_y_grid: Sign change occurs in y-grid values of LogLog\n";
      mssg +=            "interpolation.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  static const Interpolation interpolation = Interpolation::LogLog;
};

}  // namespace pndl

#endif
