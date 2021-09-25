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
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>

#include <functional>
#include <random>

namespace py = pybind11;

// LCG parameters
static std::linear_congruential_engine<uint64_t, 2806196910506780709ULL, 1,
                                       0x8000000000000000>
    lcg_engine;
static std::uniform_real_distribution<double> unit_dist;

double rang() { return unit_dist(lcg_engine); }

std::function<double()> rng(rang);

void seed(uint64_t seed) { lcg_engine.seed(seed); }

void advance_seed(unsigned long long n) { lcg_engine.discard(n); }

void init_PRNG(py::module& m) {
  m.def("seed", &seed);
  m.def("advance_seed", &advance_seed);
  m.def("rang", &rang);
  m.attr("rng") = rng;
}
