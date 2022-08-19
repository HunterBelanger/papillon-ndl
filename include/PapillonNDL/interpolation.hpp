/*
 * Papillon Nuclear Data Library
 * Copyright 2021, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Papillon Nuclear Data Library (PapillonNDL).
 *
 * PapillonNDL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PapillonNDL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PapillonNDL. If not, see <https://www.gnu.org/licenses/>.
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
#include <ostream>
#include <variant>

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

inline std::ostream& operator<<(std::ostream& out,
                                const Interpolation& interp) {
  switch (interp) {
    case Interpolation::Histogram:
      out << "Histogram";
      break;
    case Interpolation::LinLin:
      out << "LinLin";
      break;
    case Interpolation::LinLog:
      out << "LinLog";
      break;
    case Interpolation::LogLin:
      out << "LogLin";
      break;
    case Interpolation::LogLog:
      out << "LogLog";
      break;
  }

  return out;
}

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
          "Repeated values found in x-grid values of Histogram interpolation "
          "at index " +
          std::to_string(ind) + ".";
      throw PNDLException(mssg);
    }
  }

  template <class ForwardIt>
  static void verify_y_grid(ForwardIt /*first*/, ForwardIt /*last*/) {
    // No requirements
  }
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
          "Non-increasing x-grid values in LinLin interpolation.";
      throw PNDLException(mssg);
    }
  }

  template <class ForwardIt>
  static void verify_y_grid(ForwardIt /*first*/, ForwardIt /*last*/) {
    // No requirements
  }
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
          "Sign change occurs in x-grid values of LinLog interpolation.";
      throw PNDLException(mssg);
    }
  }

  template <class ForwardIt>
  static void verify_y_grid(ForwardIt /*first*/, ForwardIt /*last*/) {
    // No requirements
  }
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
          "Sign change occurs in y-grid values of LogLin interpolation.";
      throw PNDLException(mssg);
    }
  }
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
          "Sign change occurs in x-grid values of LogLog interpolation.";
      throw PNDLException(mssg);
    }
  }

  template <class ForwardIt>
  static void verify_y_grid(ForwardIt first, ForwardIt last) {
    // Make sure ther are no sign changes
    if (has_sign_change(first, last)) {
      std::string mssg =
          "Sign change occurs in y-grid values of LogLog interpolation.";
      throw PNDLException(mssg);
    }
  }
};

/**
 * @brief A generic interface for any Interpolation rule.
 */
class Interpolator {
 public:
  /**
   * @param interp Interpolation rule to use for all functions.
   */
  Interpolator(Interpolation interp) : interpolator_() {
    switch (interp) {
      case Interpolation::Histogram:
        interpolator_ = Histogram();
        break;
      case Interpolation::LinLin:
        interpolator_ = LinLin();
        break;
      case Interpolation::LinLog:
        interpolator_ = LinLog();
        break;
      case Interpolation::LogLin:
        interpolator_ = LogLin();
        break;
      case Interpolation::LogLog:
        interpolator_ = LogLog();
        break;
    }
  }

  /**
   * @brief Interpolates between (x1,y1) and (x2,y2), calculating y for a given
   * x.
   * @param x Value at which to perform the interpolation.
   * @param x1 x coordinate of the first known point.
   * @param y1 y coordinate of the first known point.
   * @param x2 x coordinate of the second known point.
   * @param y2 y coordinate of the second known point.
   */
  template <class T>
  T interpolate(T x, T x1, T y1, T x2, T y2) {
    auto doInterp = [&x, &x1, &x2, &y1, &y2](auto& interp) {
      return interp.interpolate(x, x1, y1, x2, y2);
    };
    return std::visit(doInterp, interpolator_);
  }

  /**
   * @brief Reverse interpolates between (x1,y1) and (x2,y2), calculating x for
   * a given y.
   * @param y Value at which to perform the inversion.
   * @param x1 x coordinate of the first known point.
   * @param y1 y coordinate of the first known point.
   * @param x2 x coordinate of the second known point.
   * @param y2 y coordinate of the second known point.
   */
  template <class T>
  T invert(T y, T x1, T y1, T x2, T y2) {
    auto doInvert = [&y, &x1, &x2, &y1, &y2](auto& interp) {
      return interp.invert(y, x1, y1, x2, y2);
    };
    return std::visit(doInvert, interpolator_);
  }

  /**
   * @brief Integrates between x_low and x_hi according to the provided
   *        interpolation rule, and the known points (x1,y1) and (x2,y2).
   * @param x_low Lower bound of integration.
   * @param x_hi Upper bound of integration.
   * @param x1 x coordinate of the first known point.
   * @param y1 y coordinate of the first known point.
   * @param x2 x coordinate of the second known point.
   * @param y2 y coordinate of the second known point.
   */
  template <class T>
  T integrate(T x_low, T x_hi, T x1, T y1, T x2, T y2) {
    auto doIntegrate = [&x_low, &x_hi, &x1, &x2, &y1, &y2](auto& interp) {
      return interp.integrate(x_low, x_hi, x1, y1, x2, y2);
    };
    return std::visit(doIntegrate, interpolator_);
  }

  /**
   * @brief Checks that the x-grid is valid for the given interpolation rule.
   * @param first Iterator to the first x value.
   * @param last Iterator to one after the last x value.
   */
  template <class ForwardIt>
  void verify_x_grid(ForwardIt first, ForwardIt last) {
    auto doVerifyX = [&first, &last](auto& interp) {
      interp.verify_x_grid(first, last);
    };
    std::visit(doVerifyX, interpolator_);
  }

  /**
   * @brief Checks that the y-grid is valid for the given interpolation rule.
   * @param first Iterator to the first y value.
   * @param last Iterator to one after the last y value.
   */
  template <class ForwardIt>
  void verify_y_grid(ForwardIt first, ForwardIt last) {
    auto doVerifyY = [&first, &last](auto& interp) {
      interp.verify_y_grid(first, last);
    };
    std::visit(doVerifyY, interpolator_);
  }

  /**
   * @brief Returns the current Interpolation rule.
   */
  Interpolation interpolation() const {
    if (std::holds_alternative<Histogram>(interpolator_)) {
      return Interpolation::Histogram;
    } else if (std::holds_alternative<LinLin>(interpolator_)) {
      return Interpolation::LinLin;
    } else if (std::holds_alternative<LinLog>(interpolator_)) {
      return Interpolation::LinLog;
    } else if (std::holds_alternative<LogLin>(interpolator_)) {
      return Interpolation::LogLin;
    } else {
      return Interpolation::LogLog;
    }
  }

 private:
  std::variant<Histogram, LinLin, LinLog, LogLin, LogLog> interpolator_;
};

}  // namespace pndl

#endif
