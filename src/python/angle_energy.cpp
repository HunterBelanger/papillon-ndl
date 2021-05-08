/*
 * Copyright 2021, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *
 * */
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/continuous_energy_discrete_cosines.hpp>
#include <PapillonNDL/discrete_cosines_energies.hpp>
#include <PapillonNDL/energy_angle_table.hpp>
#include <PapillonNDL/kalbach.hpp>
#include <PapillonNDL/kalbach_table.hpp>
#include <PapillonNDL/nbody.hpp>
#include <PapillonNDL/st_coherent_elastic.hpp>
#include <PapillonNDL/st_incoherent_elastic.hpp>
#include <PapillonNDL/tabular_energy_angle.hpp>
#include <PapillonNDL/uncorrelated.hpp>
#include <PapillonNDL/multiple_distribution.hpp>
#include <PapillonNDL/absorption.hpp>

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
};

void init_AngleEnergy(py::module& m) {
  py::class_<AngleEnergy, PyAngleEnergy, std::shared_ptr<AngleEnergy>>(
      m, "AngleEnergy")
      .def(py::init<>())
      .def("sample_angle_energy", &AngleEnergy::sample_angle_energy);
}

void init_Uncorrelated(py::module& m) {
  py::class_<Uncorrelated, AngleEnergy, std::shared_ptr<Uncorrelated>>(
      m, "Uncorrelated")
      .def(py::init<const AngleDistribution&,
                    std::shared_ptr<EnergyLaw>>())
      .def("sample_angle_energy", &Uncorrelated::sample_angle_energy)
      .def("angle", &Uncorrelated::angle)
      .def("energy", &Uncorrelated::energy, py::return_value_policy::reference_internal);
}

void init_NBody(py::module& m) {
  py::class_<NBody, AngleEnergy, std::shared_ptr<NBody>>(m, "NBody")
      .def(py::init<uint32_t, double, double, double>())
      .def("sample_angle_energy", &NBody::sample_angle_energy)
      .def("n", &NBody::n)
      .def("Ap", &NBody::Ap)
      .def("A", &NBody::A)
      .def("Q", &NBody::Q);
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
      .def("pdf", &KalbachTable::pdf)
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
      .def("size", &Kalbach::size);
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
      .def("pdf", &EnergyAngleTable::pdf)
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
      .def("size", &TabularEnergyAngle::size);
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
      .def("outgoing_energies", &DiscreteCosinesEnergies::outgoing_energies);
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
           &ContinuousEnergyDiscreteCosines::unit_based_interpolation);
}

void init_STCoherentElastic(py::module& m) {
  py::class_<STCoherentElastic, AngleEnergy,
             std::shared_ptr<STCoherentElastic>>(m, "STCoherentElastic")
      .def(py::init<const ACE&>())
      .def("xs", &STCoherentElastic::xs)
      .def("sample_angle_energy", &STCoherentElastic::sample_angle_energy)
      .def("bragg_edges", &STCoherentElastic::bragg_edges)
      .def("structure_factor_sum", &STCoherentElastic::structure_factor_sum);
}

void init_STInoherentElastic(py::module& m) {
  py::class_<STIncoherentElastic, AngleEnergy,
             std::shared_ptr<STIncoherentElastic>>(m, "STIncoherentElastic")
      .def(py::init<const ACE&>())
      .def("xs", py::overload_cast<>(&STIncoherentElastic::xs, py::const_), py::return_value_policy::reference_internal)
      .def("xs", py::overload_cast<double>(&STIncoherentElastic::xs, py::const_))
      .def("sample_angle_energy", &STIncoherentElastic::sample_angle_energy)
      .def("incoming_energy", &STIncoherentElastic::incoming_energy)
      .def("cosines", &STIncoherentElastic::cosines);
}

void init_MultipleDistribution(py::module& m) {
  py::class_<MultipleDistribution, AngleEnergy,
             std::shared_ptr<MultipleDistribution>>(m, "MultipleDistribution")
      .def(py::init<const std::vector<std::shared_ptr<AngleEnergy>>&,
                    const std::vector<std::shared_ptr<Tabulated1D>>&>())
      .def("sample_angle_energy", &MultipleDistribution::sample_angle_energy)
      .def("size", &MultipleDistribution::size)
      .def("distribution", &MultipleDistribution::distribution, py::return_value_policy::reference_internal)
      .def("probability", &MultipleDistribution::probability, py::return_value_policy::reference_internal);
}

void init_Absorption(py::module& m) {
  py::class_<Absorption, AngleEnergy,
             std::shared_ptr<Absorption>>(m, "Absorption")
      .def(py::init<>())
      .def("sample_angle_energy", &Absorption::sample_angle_energy);
}