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

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>

namespace pndl {

/**
 * @brief Enum to indicate the type of interpolation to use when evaluating
 *        tabulated data.
 */
enum class Interpolation : uint32_t {
  Histogram = 1, /**< y is constant in x */
  LinLin = 2,    /**< y is linear in x */
  LinLog = 3,    /**< y is linear in ln(x) */
  LogLin = 4,    /**< ln(y) is linear in x */
  LogLog = 5     /**< ln(y) is linear in ln(x) */
};

/**
 * @brief Returns true if a range of data has a sign change from positive to
 *        negative or negative to positive.
 * @param first Forward iterator to the first element in the range.
 * @param last Forward iterator to the end of the range.
 */
template <class ForwardIt>
bool has_sign_change(ForwardIt first, ForwardIt last) {
  // Get initial sign of array
  auto first_sign = std::signbit(*first);
  for (auto it = first; it != last; it++) {
    if (std::signbit(*it) != first_sign) {
      // There is a sign change in the grid !
      return true;
    }
  }
  return false;
}

/**
 * @brief Struct to perform Histogram interpolation, integration, etc.
 */
struct Histogram {
  template <class T>
  static T interpolate(T /*x*/, T /*x1*/, T y1, T /*x2*/, T /*y2*/) {
    return y1;
  }

  template <class T>
  static T invert(T /*y*/, T x1, T /*y1*/, T /*x2*/, T /*y2*/) {
    return x1;
  }

  template <class T>
  static T integrate(T x_low, T x_hi, T /*x1*/, T y1, T /*x2*/, T /*y2*/) {
    return y1 * (x_hi - x_low);
  }

  template <class ForwardIt>
  static void verify_x_grid(ForwardIt first, ForwardIt last) {
    auto it = std::adjacent_find(first, last);
    if (it != last) {
      // Grid has repeated elements which are not allowed !
      auto ind = std::distance(first, it);
      std::string mssg =
          "Histogram::verify_x_grid: Repeated values found in x-grid values\n";
      mssg +=
          "of Histogram interpolation at index " + std::to_string(ind) + ".";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  template <class ForwardIt>
  static void verify_y_grid(ForwardIt /*first*/, ForwardIt /*last*/) {
    // No requirements
  }

  static const Interpolation interpolation = Interpolation::Histogram;
};

/**
 * @brief Struct to perform LinLin interpolation, integration, etc.
 */
struct LinLin {
  template <class T>
  static T interpolate(T x, T x1, T y1, T x2, T y2) {
    return (x - x1) / (x2 - x1) * (y2 - y1) + y1;
  }

  template <class T>
  static T invert(T y, T x1, T y1, T x2, T y2) {
    return (y - y1) / (y2 - y1) * (x2 - x1) + x1;
  }

  template <class T>
  static T integrate(T x_low, T x_hi, T x1, T y1, T x2, T y2) {
    const auto numerator =
        (x_hi - x_low) * (y1 - y2) * (x_hi + x_low - 2. * x1);
    const auto denominator = 2. * (x1 - x2);
    return (numerator / denominator) + (x_hi - x_low) * y1;
  }

  template <class ForwardIt>
  static void verify_x_grid(ForwardIt first, ForwardIt last) {
    if (!std::is_sorted(first, last)) {
      // Grid elements aren't increasing, which isn't allowed for LinLin
      std::string mssg =
          "LinLin::verify_x_grid: Non-increasing x-grid values in LinLin\n";
      mssg += "interpolation.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  template <class ForwardIt>
  static void verify_y_grid(ForwardIt /*first*/, ForwardIt /*last*/) {
    // No requirements
  }

  static const Interpolation interpolation = Interpolation::LinLin;
};

/**
 * @brief Struct to perform LinLog interpolation, integration, etc.
 */
struct LinLog {
  template <class T>
  static T interpolate(T x, T x1, T y1, T x2, T y2) {
    return y1 + (y2 - y1) * std::log(x / x1) / std::log(x2 / x1);
  }

  template <class T>
  static T invert(T y, T x1, T y1, T x2, T y2) {
    return (y2 != y1) ? x1 * std::pow((x2 / x1), ((y - y1) / (y2 - y1))) : x1;
  }

  template <class T>
  static T integrate(T x_low, T x_hi, T x1, T y1, T x2, T y2) {
    const auto numerator_hi = x_hi * ((y2 - y1) * std::log(x_hi / x1) +
                                      y1 * std::log(x2 / x1) + y1 - y2);
    const auto numerator_low = x_low * ((y2 - y1) * std::log(x_low / x1) +
                                        y1 * std::log(x2 / x1) + y1 - y2);
    const auto denominator = std::log(x2 / x1);
    return (numerator_hi / denominator) - (numerator_low / denominator);
  }

  template <class ForwardIt>
  static void verify_x_grid(ForwardIt first, ForwardIt last) {
    // Make sure ther are no sign changes
    if (has_sign_change(first, last)) {
      std::string mssg =
          "LinLog::verify_x_grid: Sign change occurs in x-grid values of "
          "LinLog\n";
      mssg += "interpolation.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  template <class ForwardIt>
  static void verify_y_grid(ForwardIt /*first*/, ForwardIt /*last*/) {
    // No requirements
  }

  static const Interpolation interpolation = Interpolation::LinLog;
};

/**
 * @brief Struct to perform LogLin interpolation, integration, etc.
 */
struct LogLin {
  template <class T>
  static T interpolate(T x, T x1, T y1, T x2, T y2) {
    return y1 * std::pow((y2 / y1), (x - x1) / (x2 - x1));
  }

  template <class T>
  static T invert(T y, T x1, T y1, T x2, T y2) {
    return (y1 != y2) ? x1 + (x2 - x1) * std::log(y / y1) / std::log(y2 / y1)
                      : x1;
  }

  template <class T>
  static T integrate(T x_low, T x_hi, T x1, T y1, T x2, T y2) {
    const auto base = y2 / y1;
    const auto denominator = std::log(base);
    const auto coefficient = y1 * (x2 - x1);
    const auto exponent_hi = (x1 - x_hi) / (x1 - x2);
    const auto exponent_low = (x1 - x_low) / (x1 - x2);
    return (coefficient / denominator) *
           (std::pow(base, exponent_hi) - std::pow(base, exponent_low));
  }

  template <class ForwardIt>
  static void verify_x_grid(ForwardIt /*first*/, ForwardIt /*last*/) {
    // No requirements
  }

  template <class ForwardIt>
  static void verify_y_grid(ForwardIt first, ForwardIt last) {
    // Make sure ther are no sign changes
    if (has_sign_change(first, last)) {
      std::string mssg =
          "LogLin::verify_y_grid: Sign change occurs in y-grid values of "
          "LogLin\n";
      mssg += "interpolation.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  static const Interpolation interpolation = Interpolation::LogLin;
};

/**
 * @brief Struct to perform LogLog interpolation, integration, etc.
 */
struct LogLog {
  template <class T>
  static T interpolate(T x, T x1, T y1, T x2, T y2) {
    const auto exponent = std::log(y2 / y1) / std::log(x2 / x1);
    return y1 * std::pow(x / x1, exponent);
  }

  template <class T>
  static T invert(T y, T x1, T y1, T x2, T y2) {
    const auto exponent = std::log(y / y1) / std::log(y2 / y1);
    return x1 * std::pow(x2 / x1, exponent);
  }

  template <class T>
  static T integrate(T x_low, T x_hi, T x1, T y1, T x2, T y2) {
    const auto y2_y1 = y2 / y1;
    const auto x2_x1 = x2 / x1;
    const auto log_y2_y1 = std::log(y2_y1);
    const auto log_x2_x1 = std::log(x2_x1);
    const auto exponent = log_y2_y1 / log_x2_x1;
    const auto denominator = exponent + 1.0;
    if (std::abs(denominator) <= 1.E-12) {
      return y1 * x1 * std::log(x_hi / x_low);
    }
    return (y1 / denominator) * (x_hi * std::pow(x_hi / x1, exponent) -
                                 x_low * std::pow(x_low / x1, exponent));
  }

  template <class ForwardIt>
  static void verify_x_grid(ForwardIt first, ForwardIt last) {
    // Make sure ther are no sign changes
    if (has_sign_change(first, last)) {
      std::string mssg =
          "LogLog::verify_x_grid: Sign change occurs in x-grid values of "
          "LogLog\n";
      mssg += "interpolation.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  template <class ForwardIt>
  static void verify_y_grid(ForwardIt first, ForwardIt last) {
    // Make sure ther are no sign changes
    if (has_sign_change(first, last)) {
      std::string mssg =
          "LogLog::verify_y_grid: Sign change occurs in y-grid values of "
          "LogLog\n";
      mssg += "interpolation.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  static const Interpolation interpolation = Interpolation::LogLog;
};

}  // namespace pndl

#endif
