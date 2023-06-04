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
#include <PapillonNDL/st_coherent_elastic.hpp>
#include <PapillonNDL/st_incoherent_elastic.hpp>
#include <PapillonNDL/st_incoherent_elastic_ace.hpp>
#include <PapillonNDL/st_thermal_scattering_law.hpp>

#include "PapillonNDL/pndl_exception.hpp"

namespace pndl {

STThermalScatteringLaw::STThermalScatteringLaw(const ACE& ace,
                                               bool unit_based_interpolation)
    : zaid_(ace.zaid()),
      awr_(ace.awr()),
      temperature_(ace.temperature()),
      has_coherent_elastic_(false),
      has_incoherent_elastic_(false),
      coherent_elastic_(nullptr),
      incoherent_elastic_(nullptr),
      incoherent_inelastic_(nullptr) {
  // First try to read Incoherent Inelastic data, as all thermal scattering
  // laws have this
  try {
    incoherent_inelastic_ =
        std::make_shared<STIncoherentInelastic>(ace, unit_based_interpolation);
  } catch (PNDLException& err) {
    std::string mssg =
        "Could not construct Incoherent Inelastic scattering data.";
    err.add_to_exception(mssg);
    throw err;
  }

  // Read Coherent Elastic (if present)
  try {
    coherent_elastic_ = std::make_shared<STCoherentElastic>(ace);
  } catch (PNDLException& err) {
    std::string mssg = "Could not construct Coherent Elastic scattering data.";
    err.add_to_exception(mssg);
    throw err;
  }
  has_coherent_elastic_ = coherent_elastic_->bragg_edges().size() > 0;

  // Read Incoherent Elastic (if present)
  try {
    // Check the elastic mode, to see if this is a normal ACE file, or if
    // it is a special ACE file made by Panglos
    int32_t elastic_mode = ace.nxs(4);
    if (elastic_mode == 0) {
      incoherent_elastic_ = std::make_shared<STIncoherentElasticACE>(ace);
      has_incoherent_elastic_ = false;
    } else if (elastic_mode == 3 || elastic_mode == 5) {
      incoherent_elastic_ = std::make_shared<STIncoherentElasticACE>(ace);
      has_incoherent_elastic_ = true;
    } else if (elastic_mode == 6 && ace.jxs(6) != 0) {
      incoherent_elastic_ = std::make_shared<STIncoherentElastic>(ace);
      has_incoherent_elastic_ = true;
    } else if (elastic_mode == 6) {
      incoherent_elastic_ = std::make_shared<STIncoherentElastic>(ace);
      has_incoherent_elastic_ = false;
    } else {
      std::string mssg = "Unrecognized elastic mode in TSL ACE.";
      throw PNDLException(mssg);
    }
  } catch (PNDLException& err) {
    std::string mssg =
        "Could not construct Incoherent Elastic scattering data.";
    err.add_to_exception(mssg);
    throw err;
  }
}
}  // namespace pndl
