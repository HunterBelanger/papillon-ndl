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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/energy_grid.hpp>
#include <memory>

namespace py = pybind11;

using namespace pndl;

void init_EnergyGrid(py::module& m) {
  py::class_<EnergyGrid, std::shared_ptr<EnergyGrid>>(m, "EnergyGrid")
      .def(py::init<const ACE&, uint32_t>())
      .def(py::init<const std::vector<double>&, uint32_t>())
      .def("__getitem__", &EnergyGrid::operator[])
      .def("grid", &EnergyGrid::grid)
      .def("size", &EnergyGrid::size)
      .def("min_energy", &EnergyGrid::min_energy)
      .def("max_energy", &EnergyGrid::max_energy)
      .def("get_lower_index", &EnergyGrid::get_lower_index)
      .def("urr_min_energy", &EnergyGrid::urr_min_energy)
      .def("has_urr", &EnergyGrid::has_urr)
      .def("hash_energy_grid", &EnergyGrid::hash_energy_grid);
}
