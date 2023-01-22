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
#ifndef PAPILLON_NDL_RNG_H
#define PAPILLON_NDL_RNG_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <cstddef>
#include <cstdint>

namespace pndl {

/**
 * @brief Returns a pseudo random number on the interval [0,1).
 */
double rng();

/**
 * @brief Resets the seed of rng to the default value.
 */
void rng_reset();

/**
 * @brief Sets the seed of rng to a specific value.
 * @param seed New seed for rng.
 */
void rng_seed(std::uint64_t seed);

/**
 * @brief Advances rng by a specified number of steps.
 * @brief n Number of steps to advance rng.
 */
void rng_advance(std::uint64_t n);

}  // namespace pndl

#endif
