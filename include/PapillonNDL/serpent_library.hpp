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
#ifndef PAPILLON_SERPENT_LIBRARY_H
#define PAPILLON_SERPENT_LIBRARY_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/nd_library.hpp>
#include <string>

namespace pndl {

/**
 * @brief This class is a specialization of the NDLibrary class, which is
 *        constructed from a Serpent formated xsdir file.
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
class SerpentLibrary : public NDLibrary {
 public:
  /**
   * @param fname String with the path to the xsdir file.
   */
  SerpentLibrary(const std::string& fname);
};

}  // namespace pndl

#endif
