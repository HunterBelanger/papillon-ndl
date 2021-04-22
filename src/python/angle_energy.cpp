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
#include <PapillonNDL/tabular_energy_angle.hpp>
#include <PapillonNDL/uncorrelated.hpp>

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
      .def(py::init<std::shared_ptr<Tabulated1D>>())
      .def("sample_angle_energy", &AngleEnergy::sample_angle_energy)
      .def("probability",
           py::overload_cast<double>(&AngleEnergy::probability, py::const_))
      .def("probability",
           py::overload_cast<>(&AngleEnergy::probability, py::const_));
}

void init_Uncorrelated(py::module& m) {
  py::class_<Uncorrelated, AngleEnergy, std::shared_ptr<Uncorrelated>>(
      m, "Uncorrelated")
      .def(py::init<std::shared_ptr<AngleDistribution>,
                    std::shared_ptr<EnergyLaw>, std::shared_ptr<Tabulated1D>>())
      .def("sample_angle_energy", &Uncorrelated::sample_angle_energy)
      .def("angle", &Uncorrelated::angle)
      .def("energy", &Uncorrelated::energy);
}

void init_NBody(py::module& m) {
  py::class_<NBody, AngleEnergy, std::shared_ptr<NBody>>(m, "NBody")
      .def(py::init<uint32_t, double, double, double,
                    std::shared_ptr<Tabulated1D>>())
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
      .def(
          py::init<const std::vector<double>&, const std::vector<KalbachTable>&,
                   std::shared_ptr<Tabulated1D>>())
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
                    const std::vector<EnergyAngleTable>&,
                    std::shared_ptr<Tabulated1D>>())
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
                     &ContinuousEnergyDiscreteCosines::CEDCTable::cosines);

  py::class_<ContinuousEnergyDiscreteCosines, AngleEnergy,
             std::shared_ptr<ContinuousEnergyDiscreteCosines>>(
      m, "ContinuousEnergyDiscreteCosines")
      .def(py::init<const ACE&>())
      .def("sample_angle_energy",
           &ContinuousEnergyDiscreteCosines::sample_angle_energy)
      .def("incoming_energy", &ContinuousEnergyDiscreteCosines::incoming_energy)
      .def("size", &ContinuousEnergyDiscreteCosines::size)
      .def("tables", &ContinuousEnergyDiscreteCosines::tables)
      .def("table", &ContinuousEnergyDiscreteCosines::table);
}