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

#include <PapillonNDL/mcnp_library.hpp>
#include <PapillonNDL/nd_library.hpp>
#include <PapillonNDL/serpent_library.hpp>
#include <memory>

namespace py = pybind11;

using namespace pndl;

void init_NDLibrary(py::module& m) {
  py::class_<NDLibrary, std::shared_ptr<NDLibrary>>(m, "NDLibrary")
      .def("directory_file", &NDLibrary::directory_file)
      .def("temperatures", &NDLibrary::temperatures)
      .def("nearest_temperature", &NDLibrary::nearest_temperature)
      .def("atomic_weight_ratio", &NDLibrary::atomic_weight_ratio)
      .def("load_STNeutron", &NDLibrary::load_STNeutron, py::arg("symbol"),
           py::arg("temperature"), py::arg("tolerance") = 1.)
      .def("load_STTSL", &NDLibrary::load_STTSL, py::arg("symbol"),
           py::arg("temperature"), py::arg("tolerance") = 1.)
      .def("list_STNeutron", &NDLibrary::list_STNeutron)
      .def("list_STTSL", &NDLibrary::list_STTSL);
}

void init_MCNPLibrary(py::module& m) {
  py::class_<MCNPLibrary, NDLibrary, std::shared_ptr<MCNPLibrary>>(
      m, "MCNPLibrary")
      .def(py::init<const std::string&>());
}

void init_SerpentLibrary(py::module& m) {
  py::class_<SerpentLibrary, NDLibrary, std::shared_ptr<SerpentLibrary>>(
      m, "SerpentLibrary")
      .def(py::init<const std::string&>());
}
