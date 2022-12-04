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
#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <memory>

namespace pndl {

CrossSection::CrossSection(const ACE& ace, std::size_t i,
                           std::shared_ptr<EnergyGrid> E_grid, bool get_index,
                           bool is_heating)
    : energy_grid_(E_grid), values_(nullptr), index_(0), single_value_(false) {
  uint32_t NE = static_cast<uint32_t>(ace.nxs(2));
  if (get_index) {
    index_ = ace.xss<uint32_t>(i) - 1;
    i++;
    NE = ace.xss<uint32_t>(i);
    i++;
  }

  values_ = std::make_shared<std::vector<double>>(ace.xss(i, NE));

  if (energy_grid_->size() - index_ != values_->size()) {
    std::string mssg =
        "Different number of points in the energy grid and xs-values grid. "
        "Cross section begins at " +
        std::to_string(i) + " in XSS block.";
    throw PNDLException(mssg);
  }

  if (is_heating == false) {
    for (std::size_t l = 0; l < values_->size(); l++) {
      if ((*values_)[l] < 0.) {
        std::string mssg =
            "Negative cross section found at element " + std::to_string(l) +
            " in cross section grid starting at " + std::to_string(i) + ".";
        throw PNDLException(mssg);
      }
    }
  }
}

CrossSection::CrossSection(const std::vector<double>& xs,
                           std::shared_ptr<EnergyGrid> E_grid,
                           std::size_t index)
    : energy_grid_(E_grid),
      values_(nullptr),
      index_(static_cast<uint32_t>(index)),
      single_value_(false) {
  values_ = std::make_shared<std::vector<double>>(xs);

  if (index_ >= energy_grid_->size()) {
    std::string mssg = "Starting index is larger than size of the energy grid.";
    throw PNDLException(mssg);
  }

  for (std::size_t l = 0; l < values_->size(); l++) {
    if ((*values_)[l] < 0.) {
      std::string mssg =
          "Negative cross section found at element " + std::to_string(l) + ".";
      throw PNDLException(mssg);
    }
  }

  if (energy_grid_->size() - index_ != values_->size()) {
    std::string mssg =
        "Different number of points in the energy grid and xs-values grid.";
    throw PNDLException(mssg);
  }
}

CrossSection::CrossSection(double xs, std::shared_ptr<EnergyGrid> E_grid)
    : energy_grid_(E_grid), values_(nullptr), index_(0), single_value_(true) {
  std::vector<double> xs_tmp{xs};
  values_ = std::make_shared<std::vector<double>>(xs_tmp);

  if (values_->front() < 0.) {
    std::string mssg = "Negative cross section value provided.";
    throw PNDLException(mssg);
  }
}

std::vector<double> CrossSection::energy() const {
  return {energy_grid_->grid().begin() + static_cast<std::ptrdiff_t>(index_), energy_grid_->grid().end()};
}

}  // namespace pndl
