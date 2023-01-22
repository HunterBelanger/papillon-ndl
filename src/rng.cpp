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
#include <PapillonNDL/rng.hpp>
#include <random>

namespace pndl {

using LCG = std::linear_congruential_engine<uint64_t, 2806196910506780709ULL, 1,
                                            0x8000000000000000>;

static LCG lcg;
static std::uniform_real_distribution<double> unit_dist;

double rng() { return unit_dist(lcg); }

void rng_reset() { lcg.seed(lcg.default_seed); }

void rng_seed(std::uint64_t seed) { lcg.seed(seed); }

void rng_advance(std::uint64_t n) { lcg.discard(n); }

}  // namespace pndl
