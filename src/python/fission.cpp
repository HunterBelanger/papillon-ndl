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
#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/fission.hpp>

namespace py = pybind11;

using namespace pndl;

void init_Fission(py::module& m) {
  py::class_<Fission, std::shared_ptr<Fission>>(m, "Fission")
      .def(py::init<const ACE&, std::shared_ptr<EnergyGrid>>())
      .def(py::init<const ACE&, std::shared_ptr<EnergyGrid>, const Fission&>())
      .def("nu_total", &Fission::nu_total,
           py::return_value_policy::reference_internal)
      .def("nu_prompt", &Fission::nu_prompt,
           py::return_value_policy::reference_internal)
      .def("nu_delayed", &Fission::nu_delayed,
           py::return_value_policy::reference_internal)
      .def("prompt_spectrum", &Fission::prompt_spectrum,
           py::return_value_policy::reference_internal)
      .def("n_delayed_families", &Fission::n_delayed_families)
      .def("delayed_family", &Fission::delayed_family)
      .def("mt_list", &Fission::mt_list)
      .def("has_reaction", &Fission::has_reaction)
      .def("reaction", &Fission::reaction);
}
