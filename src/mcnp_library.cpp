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
#include <PapillonNDL/mcnp_library.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>

#include "constants.hpp"

namespace pndl {

MCNPLibrary::MCNPLibrary(const std::string& fname) : NDLibrary(fname) {
  std::filesystem::path xsdir_fname(fname);
  if (std::filesystem::exists(xsdir_fname) == false) {
    std::stringstream mssg;
    mssg << "The provided xsdir file name \"" << fname << "\" does not exist.";
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
      uint8_t Z = static_cast<uint8_t>(zaid / 1000);
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
  std::string zaid_str_new;
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
  directory_stream >> zaid_str;  // Get rid of 'directory' in stream
  bool get_zaid_str = true;
  while (directory_stream.eof() == false) {
    if (get_zaid_str) {
      directory_stream >> zaid_str;
    } else {
      zaid_str = zaid_str_new;
    }

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
      // The "ptable" keyword is on the next line. Go grab it.
      directory_stream >> ptable_str;
      get_zaid_str = true;
    } else if (ptable_str == "ptable") {
      // We just read "ptable" and it was on the main line
      get_zaid_str = true;
    } else {
      // The next item wasn't "+" or "ptable", so we actually just read the next
      // ZAID. Keep it for later.
      zaid_str_new = ptable_str;
      get_zaid_str = false;
    }

    double temp = std::stod(temp_str) * MEV_TO_EV * EV_TO_K;
    std::filesystem::path ace_path = datapath / fname_str;
    const std::regex zaid_ext_regex("([.][\\w]{3,5})");
    ACE::Type ace_type = ACE::Type::ASCII;
    const char zaid_suffix = zaid_str.back();
    if (ftype_str == "2") ace_type = ACE::Type::BINARY;
    zaid_str = std::regex_replace(zaid_str, zaid_ext_regex, "");

    if (zaid_suffix == 'c') {
      // Continuous energy neutron data
      uint32_t zaid_num = std::stoul(zaid_str);
      uint8_t Z = static_cast<uint8_t>(zaid_num / 1000);
      uint32_t A = zaid_num - (Z * 1000);
      ZAID zaid(Z, A);
      st_neutron_data_[zaid].tables.push_back({ace_path, ace_type, temp});
      st_neutron_data_[zaid].loaded_data.push_back(nullptr);
    } else if (zaid_suffix == 't') {
      // Thermal scattering law
      st_tsl_data_[zaid_str].tables.push_back({ace_path, ace_type, temp});
      st_tsl_data_[zaid_str].loaded_data.push_back(nullptr);
    }
  }

  // Entire xsdir has been read. We should now sort the entries by
  // temperature, and fill the temperature vectors.
  for (auto& lst : st_neutron_data_) {
    std::sort(lst.second.tables.begin(), lst.second.tables.end(),
              [](const TableEntry& L, const TableEntry& R) {
                return L.temperature < R.temperature;
              });

    lst.second.temperatures.resize(lst.second.tables.size(), 0.);
    for (std::size_t i = 0; i < lst.second.tables.size(); i++) {
      lst.second.temperatures[i] = lst.second.tables[i].temperature;
    }
  }

  for (auto& lst : st_tsl_data_) {
    std::sort(lst.second.tables.begin(), lst.second.tables.end(),
              [](const TableEntry& L, const TableEntry& R) {
                return L.temperature < R.temperature;
              });

    lst.second.temperatures.resize(lst.second.tables.size(), 0.);
    for (std::size_t i = 0; i < lst.second.tables.size(); i++) {
      lst.second.temperatures[i] = lst.second.tables[i].temperature;
    }
  }

  // Populate the symbol lists
  this->populate_symbol_lists();
}

}  // namespace pndl
