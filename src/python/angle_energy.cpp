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
#include "PapillonNDL/angle_energy.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/absorption.hpp>
#include <PapillonNDL/cm_distribution.hpp>
#include <PapillonNDL/continuous_energy_discrete_cosines.hpp>
#include <PapillonNDL/discrete_cosines_energies.hpp>
#include <PapillonNDL/elastic_svt.hpp>
#include <PapillonNDL/energy_angle_table.hpp>
#include <PapillonNDL/kalbach.hpp>
#include <PapillonNDL/kalbach_table.hpp>
#include <PapillonNDL/multiple_distribution.hpp>
#include <PapillonNDL/nbody.hpp>
#include <PapillonNDL/st_coherent_elastic.hpp>
#include <PapillonNDL/st_incoherent_elastic.hpp>
#include <PapillonNDL/tabular_energy_angle.hpp>
#include <PapillonNDL/uncorrelated.hpp>
#include <optional>

namespace py = pybind11;

using namespace pndl;

void init_AngleEnergyPacket(py::module& m) {
  py::class_<AngleEnergyPacket>(m, "AngleEnergyPacket")
      .def_readwrite("cosine_angle", &AngleEnergyPacket::cosine_angle)
      .def_readwrite("energy", &AngleEnergyPacket::energy);
}

// Trampoline class for abstract pndl::AngleEnergy
class PyAngleEnergy : public AngleEnergy {
 public:
  using AngleEnergy::AngleEnergy;

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override {
    PYBIND11_OVERRIDE_PURE(AngleEnergyPacket, AngleEnergy, sample_angle_energy,
                           E_in, rng);
  }

  std::optional<double> angle_pdf(double E_in, double mu) const override {
    PYBIND11_OVERRIDE_PURE(std::optional<double>, AngleEnergy, angle_pdf, E_in,
                           mu);
  }

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override {
    PYBIND11_OVERRIDE_PURE(std::optional<double>, AngleEnergy, pdf, E_in, mu,
                           E_out);
  }
};

void init_AngleEnergy(py::module& m) {
  py::class_<AngleEnergy, PyAngleEnergy, std::shared_ptr<AngleEnergy>>(
      m, "AngleEnergy")
      .def(py::init<>())
      .def("sample_angle_energy", &AngleEnergy::sample_angle_energy)
      .def("angle_pdf", &AngleEnergy::angle_pdf)
      .def("pdf", &AngleEnergy::pdf);
}

void init_Uncorrelated(py::module& m) {
  py::class_<Uncorrelated, AngleEnergy, std::shared_ptr<Uncorrelated>>(
      m, "Uncorrelated")
      .def(py::init<const AngleDistribution&, std::shared_ptr<EnergyLaw>>())
      .def("sample_angle_energy", &Uncorrelated::sample_angle_energy)
      .def("angle", &Uncorrelated::angle)
      .def("energy", &Uncorrelated::energy,
           py::return_value_policy::reference_internal)
      .def("angle_pdf", &Uncorrelated::angle_pdf)
      .def("pdf", &Uncorrelated::pdf);
}

void init_NBody(py::module& m) {
  py::class_<NBody, AngleEnergy, std::shared_ptr<NBody>>(m, "NBody")
      .def(py::init<uint32_t, double, double, double>())
      .def("sample_angle_energy", &NBody::sample_angle_energy)
      .def("n", &NBody::n)
      .def("Ap", &NBody::Ap)
      .def("A", &NBody::A)
      .def("Q", &NBody::Q)
      .def("angle_pdf", &NBody::angle_pdf)
      .def("pdf", &NBody::pdf);
}

void init_KalbachTable(py::module& m) {
  py::class_<KalbachTable>(m, "KalbachTable")
      .def(py::init<const ACE&, size_t>())
      .def(py::init<const std::vector<double>&, const std::vector<double>&,
                    const std::vector<double>&, const std::vector<double>&,
                    const std::vector<double>&, Interpolation>())
      .def("sample_energy", &KalbachTable::sample_energy)
      .def("min_energy", &KalbachTable::min_energy)
      .def("max_energy", &KalbachTable::max_energy)
      .def("energy", &KalbachTable::energy)
      .def("pdf", py::overload_cast<>(&KalbachTable::pdf, py::const_))
      .def("pdf",
           py::overload_cast<double, double>(&KalbachTable::pdf, py::const_))
      .def("angle_pdf", &KalbachTable::angle_pdf)
      .def("cdf", &KalbachTable::cdf)
      .def("R", py::overload_cast<double>(&KalbachTable::R, py::const_))
      .def("R", py::overload_cast<>(&KalbachTable::R, py::const_))
      .def("A", py::overload_cast<double>(&KalbachTable::A, py::const_))
      .def("A", py::overload_cast<>(&KalbachTable::A, py::const_))
      .def("interpolation", &KalbachTable::interpolation);
}

void init_Kalbach(py::module& m) {
  py::class_<Kalbach, AngleEnergy, std::shared_ptr<Kalbach>>(m, "Kalbach")
      .def(py::init<const std::vector<double>&,
                    const std::vector<KalbachTable>&>())
      .def("sample_angle_energy", &Kalbach::sample_angle_energy)
      .def("incoming_energy",
           py::overload_cast<>(&Kalbach::incoming_energy, py::const_))
      .def("incoming_energy",
           py::overload_cast<size_t>(&Kalbach::incoming_energy, py::const_))
      .def("table", &Kalbach::table)
      .def("size", &Kalbach::size)
      .def("angle_pdf", &Kalbach::angle_pdf)
      .def("pdf", &Kalbach::pdf);
}

void init_EnergyAngleTable(py::module& m) {
  py::class_<EnergyAngleTable>(m, "EnergyAngleTable")
      .def(py::init<const ACE&, size_t, size_t>())
      .def(py::init<const std::vector<double>&, const std::vector<double>&,
                    const std::vector<double>&, const std::vector<PCTable>&,
                    Interpolation>())
      .def(py::init<const PCTable&, const std::vector<PCTable>&>())
      .def("sample_angle_energy", &EnergyAngleTable::sample_angle_energy)
      .def("min_energy", &EnergyAngleTable::min_energy)
      .def("max_energy", &EnergyAngleTable::max_energy)
      .def("interpolation", &EnergyAngleTable::interpolation)
      .def("energy", &EnergyAngleTable::energy)
      .def("pdf", py::overload_cast<>(&EnergyAngleTable::pdf, py::const_))
      .def("pdf", py::overload_cast<double, double>(&EnergyAngleTable::pdf,
                                                    py::const_))
      .def("angle_pdf", &EnergyAngleTable::angle_pdf)
      .def("cdf", &EnergyAngleTable::cdf)
      .def("size", &EnergyAngleTable::size)
      .def("angle_table", &EnergyAngleTable::angle_table);
}

void init_TabularEnergyAngle(py::module& m) {
  py::class_<TabularEnergyAngle, AngleEnergy,
             std::shared_ptr<TabularEnergyAngle>>(m, "TabularEnergyAngle")
      .def(py::init<const std::vector<double>&,
                    const std::vector<EnergyAngleTable>&>())
      .def("sample_angle_energy", &TabularEnergyAngle::sample_angle_energy)
      .def(
          "incoming_energy",
          py::overload_cast<>(&TabularEnergyAngle::incoming_energy, py::const_))
      .def("incoming_energy",
           py::overload_cast<size_t>(&TabularEnergyAngle::incoming_energy,
                                     py::const_))
      .def("table", &TabularEnergyAngle::table)
      .def("size", &TabularEnergyAngle::size)
      .def("angle_pdf", &TabularEnergyAngle::angle_pdf)
      .def("pdf", &TabularEnergyAngle::pdf);
}

void init_DiscreteCosinesEnergies(py::module& m) {
  py::class_<DiscreteCosinesEnergies::DiscreteEnergy>(m, "DiscreteEnergy")
      .def_readwrite("energy", &DiscreteCosinesEnergies::DiscreteEnergy::energy)
      .def_readwrite("cosines",
                     &DiscreteCosinesEnergies::DiscreteEnergy::cosines);

  py::class_<DiscreteCosinesEnergies, AngleEnergy,
             std::shared_ptr<DiscreteCosinesEnergies>>(
      m, "DiscreteCosinesEnergies")
      .def(py::init<const ACE&>())
      .def("sample_angle_energy", &DiscreteCosinesEnergies::sample_angle_energy)
      .def("skewed", &DiscreteCosinesEnergies::skewed)
      .def("incoming_energy", &DiscreteCosinesEnergies::incoming_energy)
      .def("outgoing_energies", &DiscreteCosinesEnergies::outgoing_energies)
      .def("angle_pdf", &DiscreteCosinesEnergies::angle_pdf)
      .def("pdf", &DiscreteCosinesEnergies::pdf);
}

void init_ContinuousEnergyDiscreteCosines(py::module& m) {
  py::class_<ContinuousEnergyDiscreteCosines::CEDCTable>(m, "CEDCTable")
      .def_readwrite("energy",
                     &ContinuousEnergyDiscreteCosines::CEDCTable::energy)
      .def_readwrite("pdf", &ContinuousEnergyDiscreteCosines::CEDCTable::pdf)
      .def_readwrite("cdf", &ContinuousEnergyDiscreteCosines::CEDCTable::cdf)
      .def_readwrite("cosines",
                     &ContinuousEnergyDiscreteCosines::CEDCTable::cosines)
      .def("sample_energy",
           &ContinuousEnergyDiscreteCosines::CEDCTable::sample_energy);

  py::class_<ContinuousEnergyDiscreteCosines, AngleEnergy,
             std::shared_ptr<ContinuousEnergyDiscreteCosines>>(
      m, "ContinuousEnergyDiscreteCosines")
      .def(py::init<const ACE&, bool>())
      .def("sample_angle_energy",
           &ContinuousEnergyDiscreteCosines::sample_angle_energy)
      .def("incoming_energy", &ContinuousEnergyDiscreteCosines::incoming_energy)
      .def("size", &ContinuousEnergyDiscreteCosines::size)
      .def("tables", &ContinuousEnergyDiscreteCosines::tables)
      .def("table", &ContinuousEnergyDiscreteCosines::table)
      .def("unit_based_interpolation",
           &ContinuousEnergyDiscreteCosines::unit_based_interpolation)
      .def("angle_pdf", &ContinuousEnergyDiscreteCosines::angle_pdf)
      .def("pdf", &ContinuousEnergyDiscreteCosines::pdf);
}

void init_STCoherentElastic(py::module& m) {
  py::class_<STCoherentElastic, AngleEnergy,
             std::shared_ptr<STCoherentElastic>>(m, "STCoherentElastic")
      .def(py::init<const ACE&>())
      .def("xs", &STCoherentElastic::xs)
      .def("sample_angle_energy", &STCoherentElastic::sample_angle_energy)
      .def("bragg_edges", &STCoherentElastic::bragg_edges)
      .def("structure_factor_sum", &STCoherentElastic::structure_factor_sum)
      .def("angle_pdf", &STCoherentElastic::angle_pdf)
      .def("pdf", &STCoherentElastic::pdf);
}

void init_STInoherentElastic(py::module& m) {
  py::class_<STIncoherentElastic, AngleEnergy,
             std::shared_ptr<STIncoherentElastic>>(m, "STIncoherentElastic")
      .def(py::init<const ACE&>())
      .def("xs", py::overload_cast<>(&STIncoherentElastic::xs, py::const_),
           py::return_value_policy::reference_internal)
      .def("xs",
           py::overload_cast<double>(&STIncoherentElastic::xs, py::const_))
      .def("sample_angle_energy", &STIncoherentElastic::sample_angle_energy)
      .def("incoming_energy", &STIncoherentElastic::incoming_energy)
      .def("cosines", &STIncoherentElastic::cosines)
      .def("angle_pdf", &STIncoherentElastic::angle_pdf)
      .def("pdf", &STIncoherentElastic::pdf);
}

void init_MultipleDistribution(py::module& m) {
  py::class_<MultipleDistribution, AngleEnergy,
             std::shared_ptr<MultipleDistribution>>(m, "MultipleDistribution")
      .def(py::init<const std::vector<std::shared_ptr<AngleEnergy>>&,
                    const std::vector<std::shared_ptr<Tabulated1D>>&>())
      .def("sample_angle_energy", &MultipleDistribution::sample_angle_energy)
      .def("size", &MultipleDistribution::size)
      .def("distribution", &MultipleDistribution::distribution,
           py::return_value_policy::reference_internal)
      .def("probability", &MultipleDistribution::probability,
           py::return_value_policy::reference_internal)
      .def("angle_pdf", &MultipleDistribution::angle_pdf)
      .def("pdf", &MultipleDistribution::pdf);
}

void init_CMDistribution(py::module& m) {
  py::class_<CMDistribution, AngleEnergy, std::shared_ptr<CMDistribution>>(
      m, "CMDistribution")
      .def(py::init<double, double, std::shared_ptr<AngleEnergy>>())
      .def("sample_angle_energy", &CMDistribution::sample_angle_energy)
      .def("angle_pdf", &CMDistribution::angle_pdf)
      .def("pdf", &CMDistribution::pdf)
      .def("distribution", &CMDistribution::distribution,
           py::return_value_policy::reference_internal)
      .def("awr", &CMDistribution::awr)
      .def("q", &CMDistribution::q);
}

void init_Absorption(py::module& m) {
  py::class_<Absorption, AngleEnergy, std::shared_ptr<Absorption>>(m,
                                                                   "Absorption")
      .def(py::init<>())
      .def("sample_angle_energy", &Absorption::sample_angle_energy)
      .def("angle_pdf", &Absorption::angle_pdf)
      .def("pdf", &Absorption::pdf);
}

void init_ElasticSVT(py::module& m) {
  py::class_<ElasticSVT, AngleEnergy, std::shared_ptr<ElasticSVT>>(m,
                                                                   "ElasticSVT")
      .def(py::init<const AngleDistribution&, double, double, bool, double>(),
           py::arg("angle"), py::arg("awr"), py::arg("temperature"),
           py::arg("use_tar") = true, py::arg("tar_threshold") = 400.)
      .def("sample_angle_energy", &ElasticSVT::sample_angle_energy)
      .def("angle_pdf", &ElasticSVT::angle_pdf)
      .def("pdf", &ElasticSVT::pdf)
      .def("awr", &ElasticSVT::awr)
      .def("use_tar", &ElasticSVT::use_tar)
      .def("tar_threshold", &ElasticSVT::tar_threshold)
      .def("temperature", &ElasticSVT::temperature);
}
