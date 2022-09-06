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

#include <PapillonNDL/elastic.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <cmath>

#include "constants.hpp"

namespace pndl {

Elastic::Elastic(const AngleDistribution& angle, double awr, double temperature,
                 bool use_tar, double tar_threshold)
    : angle_(angle),
      awr_(awr),
      kT_(temperature * K_TO_EV * EV_TO_MEV),
      use_tar_(use_tar),
      tar_threshold_(tar_threshold) {
  if (awr_ <= 0.) {
    std::string mssg = "Atomic weight ratio must be greater than zero.";
    throw PNDLException(mssg);
  }

  if (kT_ < 0.) {
    std::string mssg = "Temperature must be greater than or equal to zero.";
    throw PNDLException(mssg);
  }

  if (tar_threshold_ < 0.) {
    std::string mssg =
        "Target At Rest threshold must be greater than or equal to zero.";
    throw PNDLException(mssg);
  }
}

double Elastic::temperature() const { return kT_ * MEV_TO_EV * EV_TO_K; }

void Elastic::set_temperature(double temperature) {
  kT_ = temperature * K_TO_EV * EV_TO_MEV;
}

void Elastic::set_tar_threshold(double tar_threshold) {
  tar_threshold_ = tar_threshold;

  if (tar_threshold_ < 0.) {
    std::string mssg =
        "Target At Rest threshold must be greater than or equal to zero.";
    throw PNDLException(mssg);
  }
}

}  // namespace pndl
