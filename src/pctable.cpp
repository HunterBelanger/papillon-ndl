/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2022, Hunter Belanger
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
#include <PapillonNDL/pctable.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>

namespace pndl {

PCTable::PCTable(const ACE& ace, std::size_t i, double normalization)
    : values_(), pdf_(), cdf_(), interp_() {
  interp_ = ace.xss<Interpolation>(i);
  if ((interp_ != Interpolation::Histogram) &&
      (interp_ != Interpolation::LinLin)) {
    std::string mssg = "Invalid interpolation of " +
                       std::to_string(static_cast<int>(interp_)) +
                       ". Index of PCTable in XSS block is " +
                       std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  uint32_t NP = ace.xss<uint32_t>(i + 1);
  if (NP == 0) {
    std::string mssg = "Cannot create a table with zero points.";
    throw PNDLException(mssg);
  }

  values_ = ace.xss(i + 2, NP);
  // Apply normalization to values
  for (auto& v : values_) v *= normalization;

  pdf_ = ace.xss(i + 2 + NP, NP);
  cdf_ = ace.xss(i + 2 + NP + NP, NP);

  if (!std::is_sorted(values_.begin(), values_.end())) {
    std::string mssg =
        "Values are not sorted. Index of PCTable in XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  if (!std::is_sorted(cdf_.begin(), cdf_.end())) {
    std::string mssg = "CDF is not sorted. Index of PCTable in XSS block is " +
                       std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  if (cdf_[cdf_.size() - 1] != 1.) {
    // If last element is close to 1, just set it to exactly 1
    if (std::abs(cdf_[cdf_.size() - 1] - 1.) < 1.E-7) {
      cdf_[cdf_.size() - 1] = 1.;
    } else {
      std::string mssg = "Last CDF entry is not 1, but " +
                         std::to_string(cdf_[cdf_.size() - 1]) +
                         ". Index of PCTable in XSS block is " +
                         std::to_string(i) + ".";
      throw PNDLException(mssg);
    }
  }

  for (const auto& p : pdf_) {
    if (p < 0.) {
      std::string mssg =
          "Negative value found in PDF. Index of PCTable in XSS block is " +
          std::to_string(i) + ".";
      throw PNDLException(mssg);
    }
  }
}

PCTable::PCTable(const std::vector<double>& values,
                 const std::vector<double>& pdf, const std::vector<double>& cdf,
                 Interpolation interp)
    : values_(values), pdf_(pdf), cdf_(cdf), interp_(interp) {
  if ((interp_ != Interpolation::Histogram) &&
      (interp_ != Interpolation::LinLin)) {
    std::string mssg = "Invalid interpolation of " +
                       std::to_string(static_cast<int>(interp_)) + ".";
    throw PNDLException(mssg);
  }

  if ((values_.size() != pdf_.size()) || (pdf_.size() != cdf_.size())) {
    std::string mssg = "Values, PDF, and CDF must have the same length.";
    throw PNDLException(mssg);
  }

  if (values_.size() == 0) {
    std::string mssg = "Cannot create a table with zero points.";
    throw PNDLException(mssg);
  }

  if (!std::is_sorted(values_.begin(), values_.end())) {
    std::string mssg = "Values are not sorted.";
    throw PNDLException(mssg);
  }

  if (!std::is_sorted(cdf_.begin(), cdf_.end())) {
    std::string mssg = "CDF is not sorted.";
    throw PNDLException(mssg);
  }

  if (cdf_[cdf_.size() - 1] != 1.) {
    // If last element is close to 1, just set it to exactly 1
    if (std::abs(cdf_[cdf_.size() - 1] - 1.) < 1.E-7) {
      cdf_[cdf_.size() - 1] = 1.;
    } else {
      std::string mssg = "Last CDF entry is not 1, but " +
                         std::to_string(cdf_[cdf_.size() - 1]) + ".";
      throw PNDLException(mssg);
    }
  }

  for (const auto& p : pdf_) {
    if (p < 0.) {
      std::string mssg = "Negative value found in PDF.";
      throw PNDLException(mssg);
    }
  }
}

}  // namespace pndl
