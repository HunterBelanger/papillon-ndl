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
#include <PapillonNDL/beta_alpha_table.hpp>
#include <PapillonNDL/pndl_exception.hpp>

namespace pndl {

BetaAlphaTable::BetaAlphaTable(const std::vector<double>& beta,
                               const std::vector<double>& pdf,
                               const std::vector<double>& cdf,
                               const std::vector<PCTable>& alpha_tables)
    : beta_(beta), pdf_(pdf), cdf_(cdf), alphas_(alpha_tables) {
  if (!std::is_sorted(beta_.begin(), beta_.end())) {
    std::string mssg = "Beta values are not sorted.";
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
