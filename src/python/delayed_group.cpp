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

#include <PapillonNDL/delayed_group.hpp>

namespace py = pybind11;

using namespace pndl;

void init_DelayedGroup(py::module& m) {
  py::class_<DelayedGroup>(m, "DelayedGroup")
      .def(py::init<const ACE&, size_t, size_t>())
      .def("decay_constant", &DelayedGroup::decay_constant)
      .def("probability", &DelayedGroup::probability,
           py::return_value_policy::reference_internal)
      .def("sample_energy", &DelayedGroup::sample_energy)
      .def("energy", &DelayedGroup::energy,
           py::return_value_policy::reference_internal);
}
