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

#include <PapillonNDL/st_incoherent_inelastic.hpp>
#include <PapillonNDL/st_thermal_scattering_law.hpp>

namespace py = pybind11;

using namespace pndl;

void init_STIncoherentInelastic(py::module& m) {
  py::class_<STIncoherentInelastic, std::shared_ptr<STIncoherentInelastic>>(
      m, "STIncoherentInelastic")
      .def(py::init<const ACE&, bool>())
      .def("cross_section", &STIncoherentInelastic::cross_section)
      .def("xs", &STIncoherentInelastic::xs)
      .def("sample_angle_energy", &STIncoherentInelastic::sample_angle_energy)
      .def("distribution", &STIncoherentInelastic::distribution)
      .def("max_energy", &STIncoherentInelastic::max_energy);
}

void init_STThermalScatteringLaw(py::module& m) {
  py::class_<STThermalScatteringLaw, std::shared_ptr<STThermalScatteringLaw>>(
      m, "STThermalScatteringLaw")
      .def(py::init<const ACE&, bool>())
      .def("zaid", &STThermalScatteringLaw::zaid)
      .def("awr", &STThermalScatteringLaw::awr)
      .def("temperature", &STThermalScatteringLaw::temperature)
      .def("max_energy", &STThermalScatteringLaw::max_energy)
      .def("has_coherent_elastic",
           &STThermalScatteringLaw::has_coherent_elastic)
      .def("has_incoherent_elastic",
           &STThermalScatteringLaw::has_incoherent_elastic)
      .def("coherent_elastic", &STThermalScatteringLaw::coherent_elastic)
      .def("incoherent_elastic", &STThermalScatteringLaw::incoherent_elastic)
      .def("incoherent_inelastic",
           &STThermalScatteringLaw::incoherent_inelastic)
      .def("coherent_elastic_xs", &STThermalScatteringLaw::coherent_elastic_xs)
      .def("incoherent_elastic_xs",
           &STThermalScatteringLaw::incoherent_elastic_xs)
      .def("incoherent_inelastic_xs",
           &STThermalScatteringLaw::incoherent_inelastic_xs);
}
