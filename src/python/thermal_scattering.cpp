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
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/st_coherent_elastic.hpp>
#include <PapillonNDL/st_incoherent_elastic_ace.hpp>
#include <PapillonNDL/st_incoherent_inelastic.hpp>
#include <PapillonNDL/st_thermal_scattering_law.hpp>
#include <PapillonNDL/st_tsl_reaction.hpp>

namespace py = pybind11;

using namespace pndl;

// Trampoline class for abstract pndl::STTSLReaction
class PySTTSLReaction : public STTSLReaction {
 public:
  using STTSLReaction::STTSLReaction;

  double xs(double E) const override {
    PYBIND11_OVERRIDE_PURE(double, STTSLReaction, xs, E);
  }

  AngleEnergyPacket sample_angle_energy(
      double E_in, const std::function<double()>& rng) const override {
    PYBIND11_OVERRIDE_PURE(AngleEnergyPacket, STTSLReaction,
                           sample_angle_energy, E_in, rng);
  }

  std::optional<double> angle_pdf(double E_in, double mu) const override {
    PYBIND11_OVERRIDE_PURE(std::optional<double>, STTSLReaction, angle_pdf,
                           E_in, mu);
  }

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override {
    PYBIND11_OVERRIDE_PURE(std::optional<double>, STTSLReaction, pdf, E_in, mu,
                           E_out);
  }
};

void init_STTSLReaction(py::module& m) {
  py::class_<STTSLReaction, PySTTSLReaction, std::shared_ptr<STTSLReaction>>(
      m, "STTSLReaction")
      .def(py::init<>())
      .def("xs", &STTSLReaction::xs);
}

void init_STCoherentElastic(py::module& m) {
  py::class_<STCoherentElastic, STTSLReaction,
             std::shared_ptr<STCoherentElastic>>(m, "STCoherentElastic")
      .def(py::init<const ACE&>())
      .def("xs", &STCoherentElastic::xs)
      .def("sample_angle_energy", &STCoherentElastic::sample_angle_energy)
      .def("bragg_edges", &STCoherentElastic::bragg_edges)
      .def("structure_factor_sum", &STCoherentElastic::structure_factor_sum)
      .def("angle_pdf", &STCoherentElastic::angle_pdf)
      .def("pdf", &STCoherentElastic::pdf);
}

void init_STInoherentElasticACE(py::module& m) {
  py::class_<STIncoherentElasticACE, AngleEnergy,
             std::shared_ptr<STIncoherentElasticACE>>(m,
                                                      "STIncoherentElasticACE")
      .def(py::init<const ACE&>())
      .def("xs", py::overload_cast<>(&STIncoherentElasticACE::xs, py::const_),
           py::return_value_policy::reference_internal)
      .def("xs",
           py::overload_cast<double>(&STIncoherentElasticACE::xs, py::const_))
      .def("sample_angle_energy", &STIncoherentElasticACE::sample_angle_energy)
      .def("incoming_energy", &STIncoherentElasticACE::incoming_energy)
      .def("cosines", &STIncoherentElasticACE::cosines)
      .def("angle_pdf", &STIncoherentElasticACE::angle_pdf)
      .def("pdf", &STIncoherentElasticACE::pdf);
}

void init_STIncoherentInelastic(py::module& m) {
  py::class_<STIncoherentInelastic, STTSLReaction,
             std::shared_ptr<STIncoherentInelastic>>(m, "STIncoherentInelastic")
      .def(py::init<const ACE&, bool>())
      .def("xs", py::overload_cast<>(&STIncoherentInelastic::xs, py::const_),
           py::return_value_policy::reference_internal)
      .def("xs",
           py::overload_cast<double>(&STIncoherentInelastic::xs, py::const_))
      .def("sample_angle_energy", &STIncoherentInelastic::sample_angle_energy)
      .def("angle_pdf", &STIncoherentInelastic::angle_pdf)
      .def("pdf", &STIncoherentInelastic::pdf)
      .def("distribution", &STIncoherentInelastic::distribution,
           py::return_value_policy::reference_internal)
      .def("max_energy", &STIncoherentInelastic::max_energy);
}

void init_STThermalScatteringLaw(py::module& m) {
  py::class_<STThermalScatteringLaw, std::shared_ptr<STThermalScatteringLaw>>(
      m, "STThermalScatteringLaw")
      .def(py::init<const ACE&, bool>(), py::arg("ace"),
           py::arg("unit_based_interpolation") = false)
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
