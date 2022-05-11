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

#include <PapillonNDL/st_incoherent_inelastic.hpp>
#include <PapillonNDL/st_thermal_scattering_law.hpp>

namespace py = pybind11;

using namespace pndl;

void init_STIncoherentInelastic(py::module& m) {
  py::class_<STIncoherentInelastic, std::shared_ptr<STIncoherentInelastic>>(
      m, "STIncoherentInelastic")
      .def(py::init<const ACE&, bool>())
      .def("xs", py::overload_cast<>(&STIncoherentInelastic::xs, py::const_),
           py::return_value_policy::reference_internal)
      .def("xs",
           py::overload_cast<double>(&STIncoherentInelastic::xs, py::const_))
      .def("sample_angle_energy", &STIncoherentInelastic::sample_angle_energy)
      .def("distribution", &STIncoherentInelastic::distribution,
           py::return_value_policy::reference_internal)
      .def("max_energy", &STIncoherentInelastic::max_energy);
}

void init_STThermalScatteringLaw(py::module& m) {
  py::class_<STThermalScatteringLaw, std::shared_ptr<STThermalScatteringLaw>>(
      m, "STThermalScatteringLaw")
      .def(py::init<const ACE&, bool>())
      .def("zaid", &STThermalScatteringLaw::zaid)
      .def("awr", &STThermalScatteringLaw::awr)
      .def("temperature", &STThermalScatteringLaw::temperature)
      .def("max_energy", &STThermalScatteringLaw::max_energy)
      .def("xs", &STThermalScatteringLaw::xs)
      .def("has_coherent_elastic",
           &STThermalScatteringLaw::has_coherent_elastic)
      .def("has_incoherent_elastic",
           &STThermalScatteringLaw::has_incoherent_elastic)
      .def("coherent_elastic", &STThermalScatteringLaw::coherent_elastic,
           py::return_value_policy::reference_internal)
      .def("incoherent_elastic", &STThermalScatteringLaw::incoherent_elastic,
           py::return_value_policy::reference_internal)
      .def("incoherent_inelastic",
           &STThermalScatteringLaw::incoherent_inelastic,
           py::return_value_policy::reference_internal);
}
