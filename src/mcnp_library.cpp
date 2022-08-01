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
#include <PapillonNDL/mcnp_library.hpp>
#include <PapillonNDL/nuclide.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>

#include "constants.hpp"

namespace pndl {

MCNPLibrary::MCNPLibrary(const std::string& fname)
    : atomic_weight_ratios_(), st_neutron_data_() {
  std::filesystem::path xsdir_fname(fname);
  if (std::filesystem::exists(xsdir_fname) == false) {
    std::stringstream mssg;
    mssg << "";
    throw PNDLException(mssg.str());
  }

  // First, read the entire file into a string
  std::stringstream xsdir_stream;
  std::ifstream xsdir_file(xsdir_fname);
  xsdir_file.seekg(0, std::ios::end);
  std::size_t xsdir_file_size = xsdir_file.tellg();
  std::string xsdir_buffer(xsdir_file_size, ' ');
  xsdir_file.seekg(0);
  xsdir_file.read(&xsdir_buffer[0], xsdir_file_size);

  //------------------------------------------------------------------
  // Now try to read the atomic weight ratios
  const std::regex awr_regex(
      "([Aa][Tt][Oo][Mm][Ii][Cc])\\s+([Ww][Ee][Ii][Gg][Hh][Tt])\\s+([Rr][Aa]["
      "Tt][Ii][Oo][Ss])(\\s+[0-9]{4,6}\\s+[0-9]+[.][0-9]+)+");
  std::smatch match;

  if (std::regex_search(xsdir_buffer, match, awr_regex) == false) {
    std::stringstream mssg;
    mssg << "Could not find atomic weight ratios in xsdir.";
    throw PNDLException(mssg.str());
  }

  // Scope reading of AWR so strings are deallocated after
  {
    std::string awr_string = match.str();
    std::stringstream awr_stream(awr_string);
    std::string input;
    uint32_t zaid;
    double awr;
    // Read three items from awr_stream, to clear out the "atomic weight ratio"
    awr_stream >> input;
    awr_stream >> input;
    awr_stream >> input;
    while (awr_stream.eof() == false) {
      awr_stream >> zaid;
      awr_stream >> awr;
      uint8_t Z = zaid / 1000;
      uint32_t A = zaid - (Z * 1000);
      ZAID zd(Z, A);
      atomic_weight_ratios_[zd] = awr;
    }
  }

  //------------------------------------------------------------------
  // Now try to read the datapath
  const std::regex datapath_regex(
      "([Dd][Aa][Tt][Aa][Pp][Aa][Tt][Hh]\\s*=?\\s*[ .\\u00C0-\\u00FF\\w/-]+)");
  std::filesystem::path datapath;
  if (std::regex_search(xsdir_buffer, match, datapath_regex) == false) {
    // We didn't find a datapath in the xsdir file, so we will assume that
    // the XS data starts from the same path as the xsdir.
    datapath = xsdir_fname.parent_path();
  } else {
    std::string datapath_string = match.str();
    const std::regex remove_datapath_regex(
        "([Dd][Aa][Tt][Aa][Pp][Aa][Tt][Hh]\\s*=?\\s*)");
    datapath_string =
        std::regex_replace(datapath_string, remove_datapath_regex, "");
    datapath = datapath_string;
  }

  //------------------------------------------------------------------
  // Now try to read the directory entries
  const std::regex directory_regex("([Dd][Ii][Rr][Ee][Cc][Tt][Oo][Rr][Yy])");
  if (std::regex_search(xsdir_buffer, match, directory_regex) == false) {
    std::stringstream mssg;
    mssg << "Could not find directory entry in xsdir.";
    throw PNDLException(mssg.str());
  }
  xsdir_buffer.erase(xsdir_buffer.begin(),
                     xsdir_buffer.begin() + match.position());
  std::stringstream directory_stream(xsdir_buffer);
  std::string zaid_str;
  std::string awr_str;
  std::string fname_str;
  std::string access_str;
  std::string ftype_str;
  std::string address_str;
  std::string tab_len_str;
  std::string record_len_str;
  std::string num_entries_str;
  std::string temp_str;
  std::string ptable_str;
  directory_stream >> zaid_str; // Get rid of 'directory' in stream
  bool get_zaid_str = true;
  while (directory_stream.eof() == false) {
    if (get_zaid_str) directory_stream >> zaid_str;

    directory_stream >> awr_str;
    directory_stream >> fname_str;

    directory_stream >> access_str;
    if (access_str == "+") {
      directory_stream >> access_str;
    }

    directory_stream >> ftype_str;
    if (ftype_str == "+") {
      directory_stream >> ftype_str;
    }

    directory_stream >> address_str;
    if (address_str == "+") {
      directory_stream >> address_str;
    }

    directory_stream >> tab_len_str;
    if (tab_len_str == "+") {
      directory_stream >> tab_len_str;
    }

    directory_stream >> record_len_str;
    if (record_len_str == "+") {
      directory_stream >> record_len_str;
    }

    directory_stream >> num_entries_str;
    if (num_entries_str == "+") {
      directory_stream >> num_entries_str;
    }

    directory_stream >> temp_str;
    if (temp_str == "+") {
      directory_stream >> temp_str;
    }

    directory_stream >> ptable_str;
    if (ptable_str == "+") {
      directory_stream >> ptable_str;
      get_zaid_str = true;
    } else if (ptable_str.size() != 6) {
      zaid_str = ptable_str;
      get_zaid_str = false;
    } else {
      get_zaid_str = true;
    }

    double temp = std::stod(temp_str) * MEV_TO_EV * EV_TO_K;
    std::filesystem::path ace_path = datapath / fname_str;
    const std::regex zaid_ext_regex("([.][\\w]{3})");
    ACE::Type ace_type = ACE::Type::ASCII;
    if (ftype_str == "2") ace_type = ACE::Type::BINARY;
    zaid_str = std::regex_replace(zaid_str, zaid_ext_regex, "");
    bool is_number = true;
    for (std::size_t i = 0; i < zaid_str.size(); i++) {
      if (std::isdigit(zaid_str[i]) == false) {
        is_number = false;
        break;
      }
    }

    if (is_number) {
      uint32_t zaid_num = std::stoul(zaid_str);
      uint8_t Z = zaid_num / 1000;
      uint32_t A = zaid_num - (Z * 1000);
      ZAID zaid(Z, A);
      st_neutron_data_[zaid].tables.push_back({ace_path, ace_type, temp});
      st_neutron_data_[zaid].loaded_data.push_back(nullptr);
    } else {
      st_tsl_data_[zaid_str].tables.push_back({ace_path, ace_type, temp});
      st_tsl_data_[zaid_str].loaded_data.push_back(nullptr);
    }
  }

  // Entire xsdir has been read. We shoudl now sort the entries by
  // temperature, and fill the temperature vectors.
  for (auto& lst : st_neutron_data_) {
    std::sort(lst.second.tables.begin(),
              lst.second.tables.end(),
              [](const TableEntry& L, const TableEntry& R) {
                return L.temperature < R.temperature;
              }
             );
    
    lst.second.temperatures.resize(lst.second.tables.size(), 0.);
    for (std::size_t i = 0; i < lst.second.tables.size(); i++) {
      lst.second.temperatures[i] = lst.second.tables[i].temperature; 
    }
  }

  for (auto& lst : st_tsl_data_) {
    std::sort(lst.second.tables.begin(),
              lst.second.tables.end(),
              [](const TableEntry& L, const TableEntry& R) {
                return L.temperature < R.temperature;
              }
             );
  
    lst.second.temperatures.resize(lst.second.tables.size(), 0.);
    for (std::size_t i = 0; i < lst.second.tables.size(); i++) {
      lst.second.temperatures[i] = lst.second.tables[i].temperature; 
    }
  }
}

const std::vector<double>& MCNPLibrary::temperatures(const std::string& symbol) const {
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

std::shared_ptr<STNeutron> MCNPLibrary::load_STNeutron(
    const std::string& symbol, double temperature, double tolerance) {
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
  for (std::size_t i = 0; i < stlist.tables.size(); i++) {
    double Tdiff = std::abs(temperature - stlist.tables[i].temperature);
    double Tdiff_1 = Tdiff + 2. * tolerance;
    if (i < stlist.tables.size() - 1) {
      Tdiff_1 = std::abs(temperature - stlist.tables[i].temperature);
    }
    if (Tdiff_1 < Tdiff) continue;

    if (Tdiff <= tolerance) {
      // We found our temperature
      // Let's check if the data is loaded
      if (stlist.loaded_data[i] == nullptr) {
        try {
          // Has yet to be loaded. We should do that. Start by loading ACE.
          ACE ace(stlist.tables[i].file, stlist.tables[i].type);

          // Now that we have the ACE, we need to construct the STNeutron
          if (stlist.first_loaded != nullptr) {
            stlist.loaded_data[i] =
                std::make_shared<STNeutron>(ace, *(stlist.first_loaded));
          } else {
            stlist.loaded_data[i] = std::make_shared<STNeutron>(ace);
            stlist.first_loaded = stlist.loaded_data[i];
          }
        } catch (PNDLException& err) {
          std::stringstream mssg;
          mssg << "Could not load STNeutron data for ACE file at ";
          mssg << stlist.tables[i].file << ".";
          err.add_to_exception(mssg.str());
          throw err;
        }
      }

      return stlist.loaded_data[i];
    }
  }

  // We didn't find a temperature within tolerance
  std::stringstream mssg;
  mssg << "Could not find data for " << symbol << " within " << tolerance;
  mssg << " Kelvin of desired temperature of " << temperature << " Kelvin.";
  throw PNDLException(mssg.str());

  // NEVER GETS HERE
  return nullptr;
}

std::shared_ptr<STThermalScatteringLaw> MCNPLibrary::load_STTSL(
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
  for (std::size_t i = 0; i < stlist.tables.size(); i++) {
    double Tdiff = std::abs(temperature - stlist.tables[i].temperature);
    double Tdiff_1 = Tdiff + 2. * tolerance;
    if (i < stlist.tables.size() - 1) {
      Tdiff_1 = std::abs(temperature - stlist.tables[i].temperature);
    }
    if (Tdiff_1 < Tdiff) continue;

    if (Tdiff <= tolerance) {
      // We found our temperature
      // Let's check if the data is loaded
      if (stlist.loaded_data[i] == nullptr) {
        try {
          // Has yet to be loaded. We should do that. Start by loading ACE.
          ACE ace(stlist.tables[i].file, stlist.tables[i].type);
          stlist.loaded_data[i] = std::make_shared<STThermalScatteringLaw>(ace);
        } catch (PNDLException& err) {
          std::stringstream mssg;
          mssg << "Could not load STThermalScatteringLaw data for ACE file at ";
          mssg << stlist.tables[i].file << ".";
          err.add_to_exception(mssg.str());
          throw err;
        }
      }

      return stlist.loaded_data[i];
    }
  }

  // We didn't find a temperature within tolerance
  std::stringstream mssg;
  mssg << "Could not find data for " << tsl_name << " within " << tolerance;
  mssg << " Kelvin of desired temperature of " << temperature << " Kelvin.";
  throw PNDLException(mssg.str());

  // NEVER GETS HERE
  return nullptr;
}

double MCNPLibrary::atomic_weight_ratio(const std::string& symbol) const {
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

ZAID MCNPLibrary::symbol_to_zaid(const std::string& symbol) const {
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
