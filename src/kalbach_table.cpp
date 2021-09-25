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
#include <PapillonNDL/kalbach_table.hpp>
#include <PapillonNDL/pndl_exception.hpp>

namespace pndl {

KalbachTable::KalbachTable(const ACE& ace, std::size_t i)
    : energy_(), pdf_(), cdf_(), R_(), A_(), interp_() {
  interp_ = ace.xss<Interpolation>(i);
  if ((interp_ != Interpolation::Histogram) &&
      (interp_ != Interpolation::LinLin)) {
    std::string mssg = "Invalid interpolation of " +
                       std::to_string(static_cast<int>(interp_)) +
                       ". Index of KalbachTable in XSS block is " +
                       std::to_string(i) + ".";
    throw PNDLException(mssg);
  }
  uint32_t NP = ace.xss<uint32_t>(i + 1);
  energy_ = ace.xss(i + 2, NP);

  pdf_ = ace.xss(i + 2 + NP, NP);
  cdf_ = ace.xss(i + 2 + NP + NP, NP);
  R_ = ace.xss(i + 2 + NP + NP + NP, NP);
  A_ = ace.xss(i + 2 + NP + NP + NP + NP, NP);

  if (!std::is_sorted(energy_.begin(), energy_.end())) {
    std::string mssg =
        "Energies are not sorted. Index of KalbachTable in XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  if (!std::is_sorted(cdf_.begin(), cdf_.end())) {
    std::string mssg =
        "CDF is not sorted. Index of KalbachTable in XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  if (cdf_[0] != 0.) {
    std::string mssg =
        "First CDF entry is not 0, but " + std::to_string(cdf_[0]) +
        ". Index of KalbachTable in XSS block is " + std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  if (cdf_[cdf_.size() - 1] != 1.) {
    // If last element is close to 1, just set it to exactly 1
    if (std::abs(cdf_[cdf_.size() - 1] - 1.) < 1.E-7) {
      cdf_[cdf_.size() - 1] = 1.;
    } else {
      std::string mssg = "Last CDF entry is not 1, but " +
                         std::to_string(cdf_[cdf_.size() - 1]) +
                         ". Index of KalbachTable in XSS block is " +
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

KalbachTable::KalbachTable(const std::vector<double>& energy,
                           const std::vector<double>& pdf,
                           const std::vector<double>& cdf,
                           const std::vector<double>& R,
                           const std::vector<double>& A, Interpolation interp)
    : energy_(energy), pdf_(pdf), cdf_(cdf), R_(R), A_(A), interp_(interp) {
  if ((interp_ != Interpolation::Histogram) &&
      (interp_ != Interpolation::LinLin)) {
    std::string mssg = "Invalid interpolation of " +
                       std::to_string(static_cast<int>(interp_)) + ".";
    throw PNDLException(mssg);
  }

  if ((energy_.size() != pdf_.size()) || (pdf_.size() != cdf_.size()) ||
      (cdf_.size() != R_.size()) || (R_.size() != A_.size())) {
    std::string mssg =
        "The outgoing energy, PDF, CDF, R, and A grids must all be the same "
        "size.";
    throw PNDLException(mssg);
  }

  if (!std::is_sorted(energy_.begin(), energy_.end())) {
    std::string mssg = "Energies are not sorted.";
    throw PNDLException(mssg);
  }

  if (!std::is_sorted(cdf_.begin(), cdf_.end())) {
    std::string mssg = "CDF is not sorted.";
    throw PNDLException(mssg);
  }

  if (cdf_[0] != 0.) {
    std::string mssg =
        "First CDF entry is not 0, but " + std::to_string(cdf_[0]) + ".";
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
