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
#include <pybind11/stl.h>

#include <PapillonNDL/reaction.hpp>
#include <PapillonNDL/reaction_base.hpp>

namespace py = pybind11;

using namespace pndl;

void init_ReactionBase(py::module& m) {
  py::class_<ReactionBase>(m, "ReactionBase")
      .def("mt", &ReactionBase::mt)
      .def("q", &ReactionBase::q)
      .def("multiplicity", &ReactionBase::yield,
           py::return_value_policy::reference_internal)
      .def("threshold", &ReactionBase::threshold)
      .def("sample_neutron_angle_energy",
           &ReactionBase::sample_neutron_angle_energy)
      .def("neutron_distribution", &ReactionBase::neutron_distribution,
           py::return_value_policy::reference_internal);
}

void init_STReaction(py::module& m) {
  py::class_<STReaction, ReactionBase>(m, "STReaction")
      .def(py::init<const ACE&, size_t, std::shared_ptr<EnergyGrid>>())
      .def(py::init<const ACE&, size_t, std::shared_ptr<EnergyGrid>,
                    const STReaction&>())
      .def(
          py::init<const CrossSection&, uint32_t, double, double, double,
                   std::shared_ptr<Function1D>, std::shared_ptr<AngleEnergy>>())
      .def("xs", &STReaction::xs);
}
