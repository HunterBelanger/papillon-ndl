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
#include <PapillonNDL/element.hpp>
#include <PapillonNDL/nd_library.hpp>
#include <PapillonNDL/nuclide.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>

#include "constants.hpp"

namespace pndl {

const std::vector<double>& NDLibrary::temperatures(
    const std::string& symbol) const {
  // first check dictionary of TSLs
  const std::regex tsl_name_regex("([\\w-]{1,6})");
  std::smatch match;
  if (std::regex_search(symbol, match, tsl_name_regex) == true) {
    std::string tsl_name = match.str();
    if (st_tsl_data_.find(tsl_name) != st_tsl_data_.end()) {
      return st_tsl_data_.at(tsl_name).temperatures;
    }
  }

  // If we didn't find a TSL, try and get a zaid
  ZAID symbol_zaid(0, 0);
  try {
    symbol_zaid = this->symbol_to_zaid(symbol);
  } catch (PNDLException& err) {
    std::stringstream mssg;
    mssg << "The symbol \"" << symbol << "\" is not a valid element or ";
    mssg << "nuclide. No thermal scattering law is associated with ";
    mssg << "this symbol.";
    err.add_to_exception(mssg.str());
    throw err;
  }

  // If we got a ZAID, check if in data map
  if (st_neutron_data_.find(symbol_zaid) == st_neutron_data_.end()) {
    // Nothing found.
    std::stringstream mssg;
    mssg << "No data associated with the symbol \"" << symbol << "\", ZAID ";
    mssg << symbol_zaid.zaid() << " was found.";
    throw PNDLException(mssg.str());
  }

  return st_neutron_data_.at(symbol_zaid).temperatures;
}

double NDLibrary::nearest_temperature(const std::string& symbol,
                                      double temperature) const {
  const std::vector<double>* temps = nullptr;

  // first check dictionary of TSLs
  const std::regex tsl_name_regex("([\\w-]{1,6})");
  std::smatch match;
  if (std::regex_search(symbol, match, tsl_name_regex) == true) {
    std::string tsl_name = match.str();
    if (st_tsl_data_.find(tsl_name) != st_tsl_data_.end()) {
      temps = &st_tsl_data_.at(tsl_name).temperatures;
    }
  }

  if (temps == nullptr) {
    // If we didn't find a TSL, try and get a zaid
    ZAID symbol_zaid(0, 0);
    try {
      symbol_zaid = this->symbol_to_zaid(symbol);
    } catch (PNDLException& err) {
      std::stringstream mssg;
      mssg << "The symbol \"" << symbol << "\" is not a valid element or ";
      mssg << "nuclide. No thermal scattering law is associated with ";
      mssg << "this symbol.";
      err.add_to_exception(mssg.str());
      throw err;
    }

    // If we got a ZAID, check if in data map
    if (st_neutron_data_.find(symbol_zaid) == st_neutron_data_.end()) {
      // Nothing found.
      std::stringstream mssg;
      mssg << "No data associated with the symbol \"" << symbol << "\", ZAID ";
      mssg << symbol_zaid.zaid() << " was found.";
      throw PNDLException(mssg.str());
    }

    temps = &st_neutron_data_.at(symbol_zaid).temperatures;
  }

  // Go through all tables to see if we find a matching temperature.
  double min_diff = temperature - (*temps)[0];
  std::size_t i_min_diff = 0;
  for (std::size_t i = 1; i < temps->size(); i++) {
    double Tdiff = std::abs(temperature - (*temps)[i]);
    if (Tdiff < min_diff) {
      min_diff = Tdiff;
      i_min_diff = i;
    }
  }

  return (*temps)[i_min_diff];
}

std::shared_ptr<STNeutron> NDLibrary::load_STNeutron(const std::string& symbol,
                                                     double temperature,
                                                     double tolerance) {
  // First get zaid of symbol
  ZAID symbol_zaid(0, 0);
  try {
    symbol_zaid = this->symbol_to_zaid(symbol);
  } catch (PNDLException& err) {
    std::stringstream mssg;
    mssg << "Could not find ZAID for provided symbol.";
    err.add_to_exception(mssg.str());
    throw err;
  }

  // If we got a ZAID, check if in data map
  if (st_neutron_data_.find(symbol_zaid) == st_neutron_data_.end()) {
    // Nothing found.
    std::stringstream mssg;
    mssg << "No data associated with the symbol \"" << symbol << "\", ZAID ";
    mssg << symbol_zaid.zaid() << " was found.";
    throw PNDLException(mssg.str());
  }

  // Get reference to the STNeutronList
  STNeutronList& stlist = st_neutron_data_[symbol_zaid];

  // Go through all tables to see if we find a matching temperature.
  double min_diff = 1.E30;
  std::size_t i_min_diff = stlist.temperatures.size();
  for (std::size_t i = 0; i < stlist.temperatures.size(); i++) {
    double Tdiff = std::abs(temperature - stlist.temperatures[i]);
    if (Tdiff <= tolerance && Tdiff < min_diff) {
      min_diff = Tdiff;
      i_min_diff = i;
    }
  }

  if (i_min_diff == stlist.temperatures.size()) {
    // We didn't find a temperature within tolerance
    std::stringstream mssg;
    mssg << "Could not find data for " << symbol << " within " << tolerance;
    mssg << " Kelvin of desired temperature of " << temperature << " Kelvin.";
    throw PNDLException(mssg.str());
  }

  // We found our temperature
  // Let's check if the data is loaded
  if (stlist.loaded_data[i_min_diff] == nullptr) {
    try {
      // Has yet to be loaded. We should do that. Start by loading ACE.
      ACE ace(stlist.tables[i_min_diff].file, stlist.tables[i_min_diff].type);

      // Now that we have the ACE, we need to construct the STNeutron
      if (stlist.first_loaded != nullptr) {
        stlist.loaded_data[i_min_diff] =
            std::make_shared<STNeutron>(ace, *(stlist.first_loaded));
      } else {
        stlist.loaded_data[i_min_diff] = std::make_shared<STNeutron>(ace);
        stlist.first_loaded = stlist.loaded_data[i_min_diff];
      }
    } catch (PNDLException& err) {
      std::stringstream mssg;
      mssg << "Could not load STNeutron data for ACE file at ";
      mssg << stlist.tables[i_min_diff].file << ".";
      err.add_to_exception(mssg.str());
      throw err;
    }
  }

  return stlist.loaded_data[i_min_diff];
}

std::shared_ptr<STThermalScatteringLaw> NDLibrary::load_STTSL(
    const std::string& symbol, double temperature, double tolerance) {
  const std::regex tsl_name_regex("([\\w-]{1,6})");
  std::smatch match;
  if (std::regex_search(symbol, match, tsl_name_regex) == false) {
    std::stringstream mssg;
    mssg << "The symbol \"" << symbol << "\" is not a valid TSL name.";
    throw PNDLException(mssg.str());
  }

  std::string tsl_name = match.str();
  if (st_tsl_data_.find(tsl_name) == st_tsl_data_.end()) {
    std::stringstream mssg;
    mssg << "Could not find \"" << tsl_name << "\" in xsdir.";
    throw PNDLException(mssg.str());
  }

  // Get reference to the STThermalScatteringLawList
  STThermalScatteringLawList& stlist = st_tsl_data_[tsl_name];

  // Go through all tables to see if we find a matching temperature.
  double min_diff = 1.E30;
  std::size_t i_min_diff = stlist.temperatures.size();
  for (std::size_t i = 0; i < stlist.temperatures.size(); i++) {
    double Tdiff = std::abs(temperature - stlist.temperatures[i]);
    if (Tdiff <= tolerance && Tdiff < min_diff) {
      min_diff = Tdiff;
      i_min_diff = i;
    }
  }

  if (i_min_diff == stlist.temperatures.size()) {
    // We didn't find a temperature within tolerance
    std::stringstream mssg;
    mssg << "Could not find data for " << tsl_name << " within " << tolerance;
    mssg << " Kelvin of desired temperature of " << temperature << " Kelvin.";
    throw PNDLException(mssg.str());
  }

  // We found our temperature
  // Let's check if the data is loaded
  if (stlist.loaded_data[i_min_diff] == nullptr) {
    try {
      // Has yet to be loaded. We should do that. Start by loading ACE.
      ACE ace(stlist.tables[i_min_diff].file, stlist.tables[i_min_diff].type);
      stlist.loaded_data[i_min_diff] =
          std::make_shared<STThermalScatteringLaw>(ace);
    } catch (PNDLException& err) {
      std::stringstream mssg;
      mssg << "Could not load STThermalScatteringLaw data for ACE file at ";
      mssg << stlist.tables[i_min_diff].file << ".";
      err.add_to_exception(mssg.str());
      throw err;
    }
  }

  return stlist.loaded_data[i_min_diff];
}

double NDLibrary::atomic_weight_ratio(const std::string& symbol) const {
  ZAID symbol_zaid(0, 0);
  try {
    symbol_zaid = this->symbol_to_zaid(symbol);
  } catch (PNDLException& err) {
    std::stringstream mssg;
    mssg << "The symbol \"" << symbol << "\" is not a valid element or ";
    mssg << "nuclide.";
    err.add_to_exception(mssg.str());
    throw err;
  }

  // If we got a ZAID, check if in data map
  if (atomic_weight_ratios_.find(symbol_zaid) == atomic_weight_ratios_.end()) {
    // Nothing found.
    std::stringstream mssg;
    mssg << "No AWR associated with the symbol \"" << symbol << "\", ZAID ";
    mssg << symbol_zaid.zaid() << " was found.";
    throw PNDLException(mssg.str());
  }

  return atomic_weight_ratios_.at(symbol_zaid);
}

ZAID NDLibrary::symbol_to_zaid(const std::string& symbol) const {
  static const std::regex is_nuclide_regex(
      "(^\\s+)?([A-Z][a-z]{0,1}[0-9]{1,3})([m][0-2])?(\\s+)?");
  static const std::regex is_element_regex("(^\\s+)?([A-Z][a-z]{0,1})(\\s+)?");

  if (std::regex_match(symbol, is_nuclide_regex)) {
    Nuclide nuc(symbol);
    return nuc.zaid();
  } else if (std::regex_match(symbol, is_element_regex)) {
    Element elem(symbol);
    return elem.zaid();
  }

  std::stringstream mssg;
  mssg << "The symbol \"" << symbol << "\" is neither a valid nuclide,";
  mssg << " nor a valid element.";
  throw PNDLException(mssg.str());

  // NEVER GETS HERE
  return ZAID(0, 0);
}

}  // namespace pndl
