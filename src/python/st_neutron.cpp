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

#include <PapillonNDL/st_neutron.hpp>
#include <memory>

namespace py = pybind11;

using namespace pndl;

void init_STNeutron(py::module& m) {
  py::class_<STNeutron, std::shared_ptr<STNeutron>>(m, "STNeutron")
      .def(py::init<const ACE&>())
      .def(py::init<const ACE&, const STNeutron&>())
      .def("zaid", &STNeutron::zaid)
      .def("awr", &STNeutron::awr)
      .def("fissile", &STNeutron::fissile)
      .def("temperature", &STNeutron::temperature)
      .def("energy_grid", &STNeutron::energy_grid)
      .def("total_xs", &STNeutron::total_xs)
      .def("elastic_xs", &STNeutron::elastic_xs)
      .def("heating_number", &STNeutron::heating_number)
      .def("fission_xs", &STNeutron::fission_xs)
      .def("disappearance_xs", &STNeutron::disappearance_xs)
      .def("photon_production_xs", &STNeutron::photon_production_xs)
      .def("mt_list", &STNeutron::mt_list)
      .def("has_reaction", &STNeutron::has_reaction)
      .def("reaction", &STNeutron::reaction)
      .def("urr_ptables", &STNeutron::urr_ptables)
      .def("elastic", py::overload_cast<>(&STNeutron::elastic, py::const_),
           py::return_value_policy::reference_internal)
      .def("elastic", py::overload_cast<>(&STNeutron::elastic),
           py::return_value_policy::reference_internal)
      .def("fission", &STNeutron::fission,
           py::return_value_policy::reference_internal)
      .def("evaluate_xs",
           py::overload_cast<double>(&STNeutron::evaluate_xs, py::const_))
      .def("evaluate_xs", py::overload_cast<double, std::size_t>(
                              &STNeutron::evaluate_xs, py::const_));
}
