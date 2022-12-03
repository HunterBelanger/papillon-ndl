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
#ifndef PANGLOS_ACE_H
#define PANGLOS_ACE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include "incoherent_inelastic.hpp"
#include "incoherent_elastic.hpp"
#include "coherent_elastic.hpp"


#include <algorithm>
#include <memory>


void write_to_ace(const LinearizedIncoherentInelastic& ii,
                  const std::unique_ptr<IncoherentElastic>& ie,
                  const std::unique_ptr<CoherentElastic>& ce,
                  double T);

#endif
