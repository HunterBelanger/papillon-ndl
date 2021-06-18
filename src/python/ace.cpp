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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/ace.hpp>

namespace py = pybind11;

using namespace pndl;

void init_ACE(py::module& m) {
  py::enum_<ACE::Type>(m, "ACEType")
      .value("ASCII", ACE::Type::ASCII)
      .value("BINARY", ACE::Type::BINARY);

  py::class_<ACE>(m, "ACE")
      .def(py::init<std::string, ACE::Type>(),
           py::arg("fname") = static_cast<std::string*>(nullptr),
           py::arg("type") = ACE::Type::ASCII)
      .def("zaid", &ACE::zaid)
      .def("temperature", &ACE::temperature)
      .def("awr", &ACE::awr)
      .def("fissile", &ACE::fissile)

      .def("izaw", py::overload_cast<size_t>(&ACE::izaw, py::const_))
      .def("izaw", py::overload_cast<size_t, size_t>(&ACE::izaw, py::const_))

      .def("nxs", py::overload_cast<size_t>(&ACE::nxs, py::const_))
      .def("nxs", py::overload_cast<size_t, size_t>(&ACE::nxs, py::const_))

      .def("jxs", py::overload_cast<size_t>(&ACE::jxs, py::const_))
      .def("jxs", py::overload_cast<size_t, size_t>(&ACE::jxs, py::const_))

      .def("xss", py::overload_cast<size_t>(&ACE::xss<double>, py::const_))
      .def("xss",
           py::overload_cast<size_t, size_t>(&ACE::xss<double>, py::const_))
      .def("zaid_id", &ACE::zaid_id)
      .def("comment", &ACE::comment)
      .def("mat", &ACE::mat)
      .def("date", &ACE::date)
      .def("save_binary", &ACE::save_binary)
      .def("ESZ", &ACE::ESZ)
      .def("NU", &ACE::NU)
      .def("MTR", &ACE::MTR)
      .def("LQR", &ACE::LQR)
      .def("TYR", &ACE::TYR)
      .def("LSIG", &ACE::LSIG)
      .def("SIG", &ACE::SIG)
      .def("LAND", &ACE::LAND)
      .def("AND", &ACE::AND)
      .def("LDLW", &ACE::LDLW)
      .def("DLW", &ACE::DLW)
      .def("DNEDL", &ACE::DNEDL)
      .def("DNED", &ACE::DNED)
      .def("DNU", &ACE::DNU)
      .def("BDD", &ACE::BDD);
}
