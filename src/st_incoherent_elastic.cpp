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
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/st_incoherent_elastic.hpp>

namespace pndl {

STIncoherentElastic::STIncoherentElastic(const ACE& ace)
    : xs_(nullptr), Nmu(0), incoming_energy_(), cosines_() {
  // Fist make sure ACE file does indeed give coherent elastic scattering
  int32_t elastic_mode = ace.nxs(4);
  if (elastic_mode == 4 || ace.jxs(3) == 0) {
    // Make a zero xs incase a user try to get the XS
    std::vector<double> E(2, 0.);
    E[1] = 100.;
    std::vector<double> xs_vals(2, 0.);
    xs_ = std::make_shared<Region1D>(E, xs_vals, Interpolation::Histogram);
  } else {
    // Get index to incident energy
    int32_t i = ace.jxs(3) - 1;
    uint32_t Ne = ace.xss<uint32_t>(i);
    incoming_energy_ = ace.xss(i + 1, Ne);
    std::vector<double> xs_vals = ace.xss(i + 1 + Ne, Ne);

    // Make sure no negative XS values
    for (size_t j = 0; j < xs_vals.size(); j++) {
      if (xs_vals[j] < 0.) {
        std::string mssg =
            "Negative cross section found at index " + std::to_string(j) + ".";
        throw PNDLException(mssg);
      }
    }

    xs_ = std::make_shared<Region1D>(incoming_energy_, xs_vals,
                                     Interpolation::LinLin);

    // Read scattering cosines
    Nmu = static_cast<uint32_t>(ace.nxs(5) + 1);
    i = ace.jxs(5) - 1;
    for (size_t ie = 0; ie < Ne; ie++) {
      std::vector<double> cosines_for_ie = ace.xss(i, Nmu);
      i += Nmu;

      // Check cosines
      if (!std::is_sorted(cosines_for_ie.begin(), cosines_for_ie.end())) {
        std::string mssg = "Cosines are not sored for incoming energy index " +
                           std::to_string(ie) + ".";
        throw PNDLException(mssg);
      }

      if (cosines_for_ie.front() < -1.) {
        std::string mssg =
            "Lowest cosine is less than -1 for incoming energy index " +
            std::to_string(ie) + ".";
        throw PNDLException(mssg);
      }

      if (cosines_for_ie.back() > 1.) {
        std::string mssg =
            "Largest cosine is greater than 1 for incoming eneergy index " +
            std::to_string(ie) + ".";
        throw PNDLException(mssg);
      }

      cosines_.push_back(cosines_for_ie);
    }  // For all incoming energies
  }
}

}  // namespace pndl
