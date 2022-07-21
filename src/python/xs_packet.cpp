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

#include <PapillonNDL/xs_packet.hpp>

namespace py = pybind11;

using namespace pndl;

void init_XSPacket(py::module& m) {
  // Initialize public member classes
  py::class_<XSPacket>(m, "XSPacket")
      .def_readwrite("total", &XSPacket::total)
      .def_readwrite("elastic", &XSPacket::elastic)
      .def_readwrite("inelastic", &XSPacket::inelastic)
      .def_readwrite("absorption", &XSPacket::absorption)
      .def_readwrite("capture", &XSPacket::capture)
      .def_readwrite("fission", &XSPacket::fission)
      .def_readwrite("heating", &XSPacket::heating)
      .def("__add__", py::overload_cast<const XSPacket&>(&XSPacket::operator+, py::const_))
      .def("__sub__", py::overload_cast<const XSPacket&>(&XSPacket::operator-, py::const_))
      .def("__mul__", &XSPacket::operator*)
      .def("__truediv__", &XSPacket::operator/)
      .def("__iadd__", &XSPacket::operator+=)
      .def("__isub__", &XSPacket::operator-=)
      .def("__imul__", &XSPacket::operator*=)
      .def("__itruediv__", &XSPacket::operator/=)
      .def("__pos__", py::overload_cast<>(&XSPacket::operator+, py::const_))
      .def("__neg__", py::overload_cast<>(&XSPacket::operator-, py::const_));
}
