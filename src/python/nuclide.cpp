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
#include <pybind11/operators.h>

#include <PapillonNDL/zaid.hpp>
#include <PapillonNDL/element.hpp>
#include <PapillonNDL/isotope.hpp>
#include <PapillonNDL/nuclide.hpp>

#include <string>

namespace py = pybind11;

using namespace pndl;

void init_ZAID(py::module& m) {
  py::class_<ZAID>(m, "ZAID")
      .def(py::init<uint8_t, uint32_t>())
      .def("Z", &ZAID::Z)
      .def("A", &ZAID::A)
      .def("zaid", &ZAID::zaid)
      .def(py::self == py::self)
      .def(py::self < py::self)
      .def("__repr__", [](const ZAID& z){ return std::to_string(z.zaid()); })
      .def("__hash__", [](const ZAID& z){ std::hash<ZAID> h; return h(z); });
}

void init_Element(py::module& m) {
  py::class_<Element>(m, "Element")
      .def(py::init<uint8_t>())
      .def("Z", &Element::Z)
      .def("atomic_numer", &Element::atomic_number)
      .def("symbol", &Element::symbol)
      .def("name", &Element::name)
      .def("zaid", &Element::zaid)
      .def(py::self == py::self)
      .def(py::self < py::self)
      .def("from_symbol", &Element::from_symbol)
      .def("from_name", &Element::from_name)
      .def("__repr__", [](const Element& e){ return e.symbol(); })
      .def("__hash__", [](const Element& e){ std::hash<Element> h; return h(e); });
}

void init_Isotope(py::module& m) {
  py::class_<Isotope>(m, "Isotope")
      .def(py::init<const Element&, uint32_t>())
      .def(py::init<uint8_t,uint32_t>())
      .def("Z", &Isotope::Z)
      .def("atomic_numer", &Isotope::atomic_number)
      .def("A", &Isotope::A)
      .def("atomic_mass", &Isotope::atomic_mass)
      .def("zaid", &Isotope::zaid)
      .def("symbol", &Isotope::symbol)
      .def("element_symbol", &Isotope::element_symbol)
      .def("element_name", &Isotope::element_name)
      .def(py::self == py::self)
      .def(py::self < py::self)
      .def("__repr__", [](const Isotope& i){ return i.symbol(); })
      .def("__hash__", [](const Isotope& i){ std::hash<Isotope> h; return h(i); });
}

void init_Nuclide(py::module& m) {
  py::class_<Nuclide>(m, "Nuclide")
      .def(py::init<const Isotope&, uint8_t>())
      .def(py::init<uint8_t,uint32_t,uint8_t>())
      .def("Z", &Nuclide::Z)
      .def("atomic_numer", &Nuclide::atomic_number)
      .def("A", &Nuclide::A)
      .def("atomic_mass", &Nuclide::atomic_mass)
      .def("level", &Nuclide::level)
      .def("zaid", &Nuclide::zaid)
      .def("symbol", &Nuclide::symbol)
      .def("isotope_symbol", &Nuclide::isotope_symbol)
      .def("element_symbol", &Nuclide::element_symbol)
      .def("element_name", &Nuclide::element_name)
      .def(py::self == py::self)
      .def(py::self < py::self)
      .def("__repr__", [](const Nuclide& n){ return n.symbol(); })
      .def("__hash__", [](const Nuclide& n){ std::hash<Nuclide> h; return h(n); });
}
