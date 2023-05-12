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
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/function_1d.hpp>
#include <PapillonNDL/linearize.hpp>

namespace py = pybind11;

using namespace pndl;

void init_Linearize(py::module& m) {
  m.def(
      "linearize",
      py::overload_cast<const std::vector<double>&, const std::vector<double>&,
                        std::function<double(double)>, double>(&linearize),
      py::arg("x"), py::arg("y"), py::arg("f"), py::arg("tolerance") = 0.001);

  m.def(
      "linearize",
      py::overload_cast<double, double, std::function<double(double)>, double>(
          &linearize),
      py::arg("x_min"), py::arg("x_max"), py::arg("f"),
      py::arg("tolerance") = 0.001);
}
