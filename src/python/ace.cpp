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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/ace.hpp>

namespace py = pybind11;

using namespace pndl;

void init_ACE(py::module& m) {
  py::enum_<ACE::Type>(m, "ACEType")
      .value("ASCII", ACE::Type::ASCII)
      .value("BINARY", ACE::Type::BINARY);

  py::class_<ACE>(m, "ACE")
      .def(py::init<std::string, ACE::Type>(), py::arg("fname"), py::arg("type") = ACE::Type::ASCII)
      .def("zaid", &ACE::zaid)
      .def("temperature", &ACE::temperature)
      .def("awr", &ACE::awr)
      .def("fissile", &ACE::fissile)

      .def("izaw", py::overload_cast<size_t>(&ACE::izaw, py::const_))
      .def("izaw", py::overload_cast<size_t, size_t>(&ACE::izaw, py::const_))

      .def("nxs", py::overload_cast<size_t>(&ACE::nxs, py::const_))
      .def("nxs", py::overload_cast<size_t, size_t>(&ACE::nxs, py::const_))

      .def("jxs", py::overload_cast<size_t>(&ACE::jxs, py::const_))
      .def("jxs", py::overload_cast<size_t, size_t>(&ACE::jxs, py::const_))

      .def("xss", py::overload_cast<size_t>(&ACE::xss<double>, py::const_))
      .def("xss",
           py::overload_cast<size_t, size_t>(&ACE::xss<double>, py::const_))
      .def("zaid_id", &ACE::zaid_id)
      .def("comment", &ACE::comment)
      .def("mat", &ACE::mat)
      .def("date", &ACE::date)
      .def("save_binary", &ACE::save_binary)
      .def("ESZ", &ACE::ESZ)
      .def("NU", &ACE::NU)
      .def("MTR", &ACE::MTR)
      .def("LQR", &ACE::LQR)
      .def("TYR", &ACE::TYR)
      .def("LSIG", &ACE::LSIG)
      .def("SIG", &ACE::SIG)
      .def("LAND", &ACE::LAND)
      .def("AND", &ACE::AND)
      .def("LDLW", &ACE::LDLW)
      .def("DLW", &ACE::DLW)
      .def("DNEDL", &ACE::DNEDL)
      .def("DNED", &ACE::DNED)
      .def("DNU", &ACE::DNU)
      .def("LUNR", &ACE::LUNR)
      .def("BDD", &ACE::BDD);
}
