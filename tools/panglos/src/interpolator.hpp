/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2023, Hunter Belanger
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
#ifndef PANGLOS_INTERPOLATOR_H
#define PANGLOS_INTERPOLATOR_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <interpolation.hpp>
using namespace njoy::interpolation;

#include <ENDFtk/types.hpp>
using namespace njoy::ENDFtk;

#include <variant>

enum class Interpolation {
  Hist = 1,
  LinLin = 2,
  LinLog = 3,
  LogLin = 4,
  LogLog = 5
};

class Interpolator {
 public:
  Interpolator(Interpolation interp) {
    switch (interp) {
      case Interpolation::Hist:
        interpr_ = Histogram();
        break;

      case Interpolation::LinLin:
        interpr_ = LinearLinear();
        break;

      case Interpolation::LinLog:
        interpr_ = LinearLogarithmic();
        break;

      case Interpolation::LogLin:
        interpr_ = LogarithmicLinear();
        break;

      case Interpolation::LogLog:
        interpr_ = LogarithmicLogarithmic();
        break;
    }
    interp_ = interp;
  }

  double interpolate(double x, double x1, double y1, double x2,
                     double y2) const {
    auto doInterp = [&x, &x1, &x2, &y1, &y2](auto& interp) {
      return interp.apply(x, x1, x2, y1, y2);
    };
    return std::visit(doInterp, interpr_);
  }

  Interpolation interpolation() const { return interp_; }

 private:
  std::variant<Histogram, LinearLinear, LinearLogarithmic, LogarithmicLinear,
               LogarithmicLogarithmic>
      interpr_;
  Interpolation interp_;
};

using Law1 =
    Table<table::Type<Histogram, table::search::Binary,
                      table::discontinuity::TakeLeft, std::vector<double>,
                      std::vector<double>>,
          table::left::interval::Throws, table::right::interval::Throws>;

using Law2 =
    Table<table::Type<LinearLinear, table::search::Binary,
                      table::discontinuity::TakeLeft, std::vector<double>,
                      std::vector<double>>,
          table::left::interval::Throws, table::right::interval::Throws>;

using Law3 =
    Table<table::Type<LinearLogarithmic, table::search::Binary,
                      table::discontinuity::TakeLeft, std::vector<double>,
                      std::vector<double>>,
          table::left::interval::Throws, table::right::interval::Throws>;

using Law4 =
    Table<table::Type<LogarithmicLinear, table::search::Binary,
                      table::discontinuity::TakeLeft, std::vector<double>,
                      std::vector<double>>,
          table::left::interval::Throws, table::right::interval::Throws>;

using Law5 =
    Table<table::Type<LogarithmicLogarithmic, table::search::Binary,
                      table::discontinuity::TakeLeft, std::vector<double>,
                      std::vector<double>>,
          table::left::interval::Throws, table::right::interval::Throws>;

using ENDFvariant = Table<table::Variant<Law1, Law2, Law3, Law4, Law5>>;
using Tab1 = Table<table::Vector<ENDFvariant>>;

// Function to generate a Tab1 from the breakpoints, interpolations,
// x array, and y array
inline Tab1 makeTab1(AllRange<long> breakpoints, AllRange<long> interpolations,
                     AllRange<double> x, AllRange<double> y) {
  std::vector<ENDFvariant> core;

  auto brkSize = ranges::size(breakpoints);
  size_t low = 0;
  size_t hi = 0;

  for (size_t i = 0; i < brkSize; i++) {
    int law = *(interpolations.begin() + i);
    hi = *(breakpoints.begin() + i);

    std::vector<double> xInterval = {x.begin() + low, x.begin() + hi};
    std::vector<double> yInterval = {y.begin() + low, y.begin() + hi};

    if (law == 1) {
      core.push_back(Law1(std::move(xInterval), std::move(yInterval)));
    } else if (law == 2) {
      core.push_back(Law2(std::move(xInterval), std::move(yInterval)));
    } else if (law == 3) {
      core.push_back(Law3(std::move(xInterval), std::move(yInterval)));
    } else if (law == 4) {
      core.push_back(Law4(std::move(xInterval), std::move(yInterval)));
    } else if (law == 5) {
      core.push_back(Law5(std::move(xInterval), std::move(yInterval)));
    } else {
      // TODO unknown interpolation law
      throw std::runtime_error("Unkown interpolation law " +
                               std::to_string(law));
    }

    low = hi - 1;

    // Check for discontinuity
    if (low < ranges::size(x) - 1) {
      if (*(x.begin() + low) == *(x.begin() + low + 1)) {
        low++;
      }
    }
  }

  Tab1 table(std::move(core));

  return table;
}

#endif
