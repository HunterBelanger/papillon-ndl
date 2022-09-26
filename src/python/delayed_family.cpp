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
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>

#include <PapillonNDL/delayed_family.hpp>

namespace py = pybind11;

using namespace pndl;

void init_DelayedFamily(py::module& m) {
  py::class_<DelayedFamily>(m, "DelayedFamily")
      .def(py::init<const ACE&, size_t, size_t>())
      .def("decay_constant", &DelayedFamily::decay_constant)
      .def("probability", &DelayedFamily::probability,
           py::return_value_policy::reference_internal)
      .def("sample_energy", &DelayedFamily::sample_energy)
      .def("energy", &DelayedFamily::energy,
           py::return_value_policy::reference_internal);
}
