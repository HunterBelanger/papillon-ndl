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
#include <pybind11/stl.h>

#include <PapillonNDL/angle_table.hpp>
#include <PapillonNDL/equiprobable_angle_bins.hpp>
#include <PapillonNDL/isotropic.hpp>
#include <PapillonNDL/legendre.hpp>

namespace py = pybind11;

using namespace pndl;

// Trampoline class for abstract pndl::AngleLaw
class PyAngleLaw : public AngleLaw {
 public:
  using AngleLaw::AngleLaw;

  double sample_mu(const std::function<double()>& rng) const override {
    PYBIND11_OVERRIDE_PURE(double, AngleLaw, sample_mu, rng);
  }

  double pdf(double mu) const override {
    PYBIND11_OVERRIDE_PURE(double, AngleLaw, pdf, mu);
  }
};

void init_AngleLaw(py::module& m) {
  py::class_<AngleLaw, PyAngleLaw, std::shared_ptr<AngleLaw>>(m, "AngleLaw")
      .def(py::init<>())
      .def("sample_mu", &AngleLaw::sample_mu)
      .def("pdf", &AngleLaw::pdf);
}

void init_Isotropic(py::module& m) {
  py::class_<Isotropic, AngleLaw, std::shared_ptr<Isotropic>>(m, "Isotropic")
      .def(py::init<>())
      .def("sample_mu", &Isotropic::sample_mu)
      .def("pdf", &Isotropic::pdf);
}

void init_EquiprobableAngleBins(py::module& m) {
  py::class_<EquiprobableAngleBins, AngleLaw,
             std::shared_ptr<EquiprobableAngleBins>>(m, "EquiprobableAngleBins")
      .def(py::init<const ACE&, size_t>())
      .def(py::init<const std::vector<double>&>())
      .def("sample_mu", &EquiprobableAngleBins::sample_mu)
      .def("pdf", &EquiprobableAngleBins::pdf)
      .def("size", &EquiprobableAngleBins::size)
      .def("bin_bounds", &EquiprobableAngleBins::bin_bounds);
}

void init_AngleTable(py::module& m) {
  py::class_<AngleTable, AngleLaw, std::shared_ptr<AngleTable>>(m, "AngleTable")
      .def(py::init<const ACE&, size_t>())
      .def(py::init<const std::vector<double>&, const std::vector<double>&,
                    const std::vector<double>&, Interpolation>())
      .def(py::init<const Legendre&>())
      .def(py::init<const PCTable&>())
      .def("sample_mu", &AngleTable::sample_mu)
      .def("size", &AngleTable::size)
      .def("cosines", &AngleTable::cosines)
      .def("pdf", py::overload_cast<>(&AngleTable::pdf, py::const_))
      .def("pdf", py::overload_cast<double>(&AngleTable::pdf, py::const_))
      .def("cdf", &AngleTable::cdf)
      .def("interpolation", &AngleTable::interpolation);
}

void init_Legendre(py::module& m) {
  py::class_<Legendre, AngleLaw, std::shared_ptr<Legendre>>(m, "Legendre")
      .def(py::init<>())
      .def(py::init<const std::vector<double>&>())
      .def("sample_mu", &Legendre::sample_mu)
      .def("pdf", &Legendre::pdf)
      .def("set_moment", &Legendre::set_moment)
      .def("coefficients", &Legendre::coefficients);
}
