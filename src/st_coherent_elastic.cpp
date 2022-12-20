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
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/st_coherent_elastic.hpp>

#include "constants.hpp"

namespace pndl {

STCoherentElastic::STCoherentElastic(const ACE& ace)
    : bragg_edges_(), structure_factor_sum_() {
  // Fist make sure ACE file does indeed give coherent elastic scattering
  int32_t elastic_mode = ace.nxs(4);
  if (elastic_mode == 4 || elastic_mode == 5) {
    // Get index to Bragg edge and structure data
    std::size_t i = static_cast<std::size_t>(ace.jxs(3) - 1);
    uint32_t Ne = ace.xss<uint32_t>(i);
    bragg_edges_ = ace.xss(i + 1, Ne);
    structure_factor_sum_ = ace.xss(i + 1 + Ne, Ne);

    // Make sure Bragg edges are all positive and sorted
    if (!std::is_sorted(bragg_edges_.begin(), bragg_edges_.end())) {
      std::string mssg = "Bragg edges are not sorted.";
      throw PNDLException(mssg);
    }

    if (bragg_edges_.front() < 0.) {
      std::string mssg = "Negative Bragg edges found.";
      throw PNDLException(mssg);
    }
  }
}

}  // namespace pndl
