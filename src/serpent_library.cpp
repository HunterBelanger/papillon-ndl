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
#include <PapillonNDL/serpent_library.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_set>

namespace pndl {

SerpentLibrary::SerpentLibrary(const std::string& fname)
    : NDLibrary() {
  std::filesystem::path xsdir_fname(fname);
  if (std::filesystem::exists(xsdir_fname) == false) {
    std::stringstream mssg;
    mssg << "The provided xsdir file name \"" << fname << "\" ";
    mssg << "does not exist.";
    throw PNDLException(mssg.str());
  }

  // First, read the entire file into a string
  std::ifstream xsdir_file(xsdir_fname);
  xsdir_file.seekg(0, std::ios::end);
  std::size_t xsdir_file_size = xsdir_file.tellg();
  std::string xsdir_buffer(xsdir_file_size, ' ');
  xsdir_file.seekg(0);
  xsdir_file.read(&xsdir_buffer[0], xsdir_file_size);
  std::stringstream xsdir_stream(xsdir_buffer);

  // Read the xsdir line by line
  std::string alias_str;
  std::string zaid_str;
  std::string type_str;
  std::string ZA_str;
  std::string I_str;
  std::string AW_str;
  std::string T_str;
  std::string bin_str;
  std::string path_str;
  std::unordered_set<std::string> read_zaids;
  while (xsdir_stream.eof() == false) {
    xsdir_stream >> alias_str;
    xsdir_stream >> zaid_str;
    xsdir_stream >> type_str;
    xsdir_stream >> ZA_str;
    xsdir_stream >> I_str;
    xsdir_stream >> AW_str;
    xsdir_stream >> T_str;
    xsdir_stream >> bin_str;
    xsdir_stream >> path_str; 

    // If this isn't a neutron or tsl file, read next line
    if (type_str[0] != '1' && type_str[0] != '3') continue;

    // Check if this is an alias. If so, read next line
    if (read_zaids.find(zaid_str) != read_zaids.end()) continue;
    read_zaids.emplace(zaid_str);
    
    ACE::Type ace_type = ACE::Type::ASCII;
    if (bin_str[0] == '1') ace_type = ACE::Type::BINARY;
    double temp = std::stod(T_str);
    std::filesystem::path ace_path = path_str;
    const std::regex zaid_ext_regex("([.][\\w]{3})");
    zaid_str = std::regex_replace(zaid_str, zaid_ext_regex, "");
    bool is_number = true;
    for (std::size_t i = 0; i < zaid_str.size(); i++) {
      if (std::isdigit(zaid_str[i]) == false) {
        is_number = false;
        break;
      }
    }
    
    if (is_number) {
      // Free-gass Neutron data
      uint32_t ZA = std::stoul(ZA_str);
      uint8_t Z = ZA / 1000;
      uint32_t A = ZA - (Z * 1000);
      ZAID zaid(Z, A);
      st_neutron_data_[zaid].tables.push_back({ace_path, ace_type, temp});
      st_neutron_data_[zaid].loaded_data.push_back(nullptr);
      double awr = std::stod(AW_str);
      atomic_weight_ratios_[zaid] = awr;
    } else {
      // Thermal Scattering Law 
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

}  // namespace pndl
