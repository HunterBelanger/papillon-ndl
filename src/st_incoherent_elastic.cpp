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
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/st_incoherent_elastic.hpp>
#include <string>

namespace pndl {

STIncoherentElastic::STIncoherentElastic(const ACE& ace) : xs_(-1.), W_(0.) {
  const int32_t elastic_mode = ace.nxs(4);
  if (elastic_mode == 3 || elastic_mode == 5) {
    std::string mssg = "Cannot construct from ACE TSL file.";
    throw PNDLException(mssg);
  } else if (elastic_mode == 6 && ace.jxs(6) != 0) {
    std::size_t i = static_cast<std::size_t>(ace.jxs(6) - 1);
    xs_ = ace.xss(i);
    W_ = ace.xss(i + 1);

    if (xs_ <= 0.) {
      std::string mssg = "Characteristic bound cross section is <= 0.";
      throw PNDLException(mssg);
    }

    if (W_ <= 0.) {
      std::string mssg = "Debye-Waller integral is <= 0.";
      throw PNDLException(mssg);
    }
  }
}

}  // namespace pndl
