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

#include <PapillonNDL/ce_neutron.hpp>
#include <PapillonNDL/ce_neutron_base.hpp>

namespace py = pybind11;

using namespace pndl;

void init_CENeutronBase(py::module& m) {
  py::class_<CENeutronBase>(m, "CENeutronBase")
      .def("zaid", &CENeutronBase::zaid)
      .def("awr", &CENeutronBase::awr)
      .def("fissile", &CENeutronBase::fissile)
      .def("elastic_angle_distribution",
           &CENeutronBase::elastic_angle_distribution)
      .def("nu_total", &CENeutronBase::nu_total,
           py::return_value_policy::reference_internal)
      .def("nu_prompt", &CENeutronBase::nu_prompt,
           py::return_value_policy::reference_internal)
      .def("nu_delayed", &CENeutronBase::nu_delayed,
           py::return_value_policy::reference_internal)
      .def("n_delayed_groups", &CENeutronBase::n_delayed_groups)
      .def("delayed_group", &CENeutronBase::delayed_group)
      .def("mt_list", &CENeutronBase::mt_list)
      .def("has_reaction", &CENeutronBase::has_reaction);
}

void init_STNeutron(py::module& m) {
  py::class_<STNeutron, CENeutronBase>(m, "STNeutron")
      .def(py::init<const ACE&>())
      .def(py::init<const ACE&, const STNeutron&>())
      .def("temperature", &STNeutron::temperature)
      .def("energy_grid", &STNeutron::energy_grid)
      .def("total_xs", &STNeutron::total_xs)
      .def("elastic_xs", &STNeutron::elastic_xs)
      .def("fission_xs", &STNeutron::fission_xs)
      .def("disappearance_xs", &STNeutron::disappearance_xs)
      .def("photon_production_xs", &STNeutron::photon_production_xs)
      .def("reaction", &STNeutron::reaction);
}
