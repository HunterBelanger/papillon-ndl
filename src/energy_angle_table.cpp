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
#include <PapillonNDL/energy_angle_table.hpp>
#include <PapillonNDL/pndl_exception.hpp>

namespace pndl {

EnergyAngleTable::EnergyAngleTable(const ACE& ace, std::size_t i,
                                   std::size_t JED)
    : energy_(), pdf_(), cdf_(), angles_(), interp_() {
  interp_ = ace.xss<Interpolation>(i);
  if ((interp_ != Interpolation::Histogram) &&
      (interp_ != Interpolation::LinLin)) {
    std::string mssg = "Invalid interpolation of " +
                       std::to_string(static_cast<int>(interp_)) +
                       ". Index of EnergyAngleTable in XSS block is " +
                       std::to_string(i) + ".";
    throw PNDLException(mssg);
  }
  uint32_t NP = ace.xss<uint32_t>(i + 1);
  energy_ = ace.xss(i + 2, NP);

  pdf_ = ace.xss(i + 2 + NP, NP);
  cdf_ = ace.xss(i + 2 + NP + NP, NP);

  if (!std::is_sorted(energy_.begin(), energy_.end())) {
    std::string mssg =
        "Energies are not sorted. Index of EnergyAngleTable in XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  if (!std::is_sorted(cdf_.begin(), cdf_.end())) {
    std::string mssg =
        "CDF is not sorted. Index of EnergyAngleTable in XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  if (cdf_[0] != 0.) {
    std::string mssg = "First CDF entry is not 0, but " +
                       std::to_string(cdf_[0]) +
                       ". Index of EnergyAngleTable in XSS block is " +
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
                         ". Index of EnergyAngleTable in XSS block is " +
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

  std::vector<int32_t> locs = ace.xss<int32_t>(i + 2 + NP + NP + NP, NP);
  for (std::size_t j = 0; j < locs.size(); j++) {
    int32_t loc = locs[j];
    std::size_t l = JED + std::abs(loc) - 1;
    try {
      angles_.emplace_back(ace, l);
    } catch (PNDLException& error) {
      std::string mssg = "Couldn't create angle table for " +
                         std::to_string(j) + "th energy " +
                         std::to_string(energy_[j]) +
                         " MeV. Index of EnergyAngleTable in XSS block is " +
                         std::to_string(i) + ".";
      error.add_to_exception(mssg);
      throw error;
    }
  }
}

EnergyAngleTable::EnergyAngleTable(const std::vector<double>& outgoing_energy,
                                   const std::vector<double>& pdf,
                                   const std::vector<double>& cdf,
                                   const std::vector<PCTable>& angle_tables,
                                   Interpolation interp)
    : energy_(outgoing_energy),
      pdf_(pdf),
      cdf_(cdf),
      angles_(angle_tables),
      interp_(interp) {
  if ((interp_ != Interpolation::Histogram) &&
      (interp_ != Interpolation::LinLin)) {
    std::string mssg = "Invalid interpolation of " +
                       std::to_string(static_cast<int>(interp_)) + ".";
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

EnergyAngleTable::EnergyAngleTable(const PCTable& outgoing_energy,
                                   const std::vector<PCTable>& angle_tables)
    : energy_(outgoing_energy.values()),
      pdf_(outgoing_energy.pdf()),
      cdf_(outgoing_energy.cdf()),
      angles_(angle_tables),
      interp_(outgoing_energy.interpolation()) {}

}  // namespace pndl
