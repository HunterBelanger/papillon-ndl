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
#ifndef PAPILLON_ND_LIBRARY_H
#define PAPILLON_ND_LIBRARY_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/st_neutron.hpp>
#include <PapillonNDL/st_thermal_scattering_law.hpp>
#include <PapillonNDL/zaid.hpp>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace pndl {

/**
 * @brief This class acts as an interface to easily retrieve nuclear data
 *        for the desried isotope or thermal scattering law, at the desired
 *        temperature. This interface makes sure that data is only ever loaded
 *        from an ACE file once, and all scattering distributions for all
 *        STNeutron instances of the same nuclide are shared, conserving memory.
 *
 * @warning Due to historical reasons, in many ACE libraries (mostly those
 *          distributed for use with MCNP by LANL), Am242m1 has been given
 *          a ZAID of 95242, and Am242 has been given a ZAID of 95642. If
 *          this is the case for the ACE library you are using, then these
 *          two evaluations will be switched, and looking up the symbol
 *          "Am242" will actually provide the evaluation for Am242m1. This
 *          can be corrected by modifying the xsdir file itself, swaping
 *          the ZAID identifiers for the two evaluations.
 */
class NDLibrary {
 public:
  virtual ~NDLibrary() = default;

  /**
   * @brief Returns the atomic weight ratio for the nuclide associated with
   *        symbol, which was tabulated in the xsdir file.
   * @param symbol String contianing the symbol for the nucldie.
   */
  double atomic_weight_ratio(const std::string& symbol) const;

  /**
   * @breif Returns a reference to a vector of all temperatures for which
   *        data is given for the provided symbol.
   * @param symbol String containing the symbol for the nuclide, or the name
   *               of the thermal scattering law.
   */
  const std::vector<double>& temperatures(const std::string& symbol) const;

  /**
   * @brief Returns the closest nearest provided temperature for the given
   *        symbol-temperature combination.
   * @param symbol String contianing the symbol for the nuclide, or the name
   *               of the thermal scattering law to be queried.
   * @param temperature Desired temperature for data.
   */
  double nearest_temperature(const std::string& symbol,
                             double temperature) const;

  /**
   * @brief Returns a shared pointer to an STNeutron instance.
   * @param symbol String for the symbol of the desired nuclide.
   * @param temperature Desired temperature for the returned STNeutron data
   *                    in Kelvin.
   * @param tolerance Temperature tolerance in Kelvin for the returned data.
   *                  If no data is found at the desired temperature, the data
   *                  at the closest temerature, within +/- tolerance will be
   *                  returned instead.
   */
  std::shared_ptr<STNeutron> load_STNeutron(const std::string& symbol,
                                            double temperature,
                                            double tolerance = 1.);

  /**
   * @brief Returns a shared pointer to an STThermalScatteringLaw instance.
   * @param symbol String for the name of the desired thermal scattering law.
   * @param temperature Desired temperature for the returned scattering law
   *                    in Kelvin.
   * @param tolerance Temperature tolerance in Kelvin for the returned data.
   *                  If no data is found at the desired temperature, the data
   *                  at the closest temerature, within +/- tolerance will be
   *                  returned instead.
   */
  std::shared_ptr<STThermalScatteringLaw> load_STTSL(const std::string& symbol,
                                                     double temperature,
                                                     double tolerance = 1.);

  /**
   * @breif Returns a vector containing all of the available symbols for
   *        STNeutron data.
   */
  const std::vector<std::string>& list_STNeutron() const {
    return st_neutron_symbols_;
  }

  /**
   * @breif Returns a vector containing all of the available symbols for
   *        STThermalScatteringLaw data.
   */
  const std::vector<std::string>& list_STTSL() const { return st_tsl_symbols_; }

  /**
   * @brief Returns the path to the file containing the library directory.
   */
  const std::string& directory_file() const { return xsdir_fname_; }

 protected:
  NDLibrary(const std::string& fname)
      : xsdir_fname_(fname),
        atomic_weight_ratios_(),
        st_neutron_data_(),
        st_tsl_data_(),
        st_neutron_symbols_(),
        st_tsl_symbols_(){};

  struct TableEntry {
    std::filesystem::path file;
    ACE::Type type;
    double temperature;
  };

  struct STNeutronList {
    std::vector<TableEntry> tables;
    std::vector<std::shared_ptr<STNeutron>> loaded_data;
    std::vector<double> temperatures;
    std::shared_ptr<STNeutron> first_loaded;
  };

  struct STThermalScatteringLawList {
    std::vector<TableEntry> tables;
    std::vector<std::shared_ptr<STThermalScatteringLaw>> loaded_data;
    std::vector<double> temperatures;
  };

  std::string xsdir_fname_;
  std::unordered_map<ZAID, double> atomic_weight_ratios_;
  std::unordered_map<ZAID, STNeutronList> st_neutron_data_;
  std::unordered_map<std::string, STThermalScatteringLawList> st_tsl_data_;
  std::vector<std::string> st_neutron_symbols_;
  std::vector<std::string> st_tsl_symbols_;

  ZAID symbol_to_zaid(const std::string& symbol) const;
  void populate_symbol_lists();
};

}  // namespace pndl

#endif
