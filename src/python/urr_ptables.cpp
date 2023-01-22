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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/urr_ptables.hpp>
#include <memory>

namespace py = pybind11;

using namespace pndl;

void init_URRPtable(py::module& m) {
  // Initialize public member classes
  py::class_<URRPTables::PTable::XSBand>(m, "XSBand")
      .def_readwrite("total", &URRPTables::PTable::XSBand::total)
      .def_readwrite("elastic", &URRPTables::PTable::XSBand::elastic)
      .def_readwrite("fission", &URRPTables::PTable::XSBand::fission)
      .def_readwrite("capture", &URRPTables::PTable::XSBand::capture)
      .def_readwrite("heating", &URRPTables::PTable::XSBand::heating);

  py::class_<URRPTables::PTable>(m, "PTable")
      .def_readwrite("cdf", &URRPTables::PTable::cdf)
      .def_readwrite("xs_bands", &URRPTables::PTable::xs_bands);

  py::class_<URRPTables, std::shared_ptr<URRPTables>>(m, "URRPTables")
      .def(py::init<const ACE&, const CrossSection&, const CrossSection&,
                    const CrossSection&, const CrossSection&,
                    const CrossSection&, const CrossSection&,
                    const std::vector<STReaction>&>())
      .def("is_valid", &URRPTables::is_valid)
      .def("evaluate_xs", py::overload_cast<double, std::size_t, double>(
                              &URRPTables::evaluate_xs, py::const_))
      .def("evaluate_xs", py::overload_cast<double, double>(
                              &URRPTables::evaluate_xs, py::const_))
      .def("min_energy", &URRPTables::min_energy)
      .def("max_energy", &URRPTables::max_energy)
      .def("energy_in_range", &URRPTables::energy_in_range)
      .def("energy", &URRPTables::energy)
      .def("ptables", &URRPTables::ptables)
      .def("n_xs_bands", &URRPTables::n_xs_bands)
      .def("xs_factors", &URRPTables::xs_factors)
      .def("inelastic_competition", &URRPTables::inelastic_competition)
      .def("absorption_competition", &URRPTables::absorption_competition);
}
