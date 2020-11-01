/*
 * Copyright 2020, Hunter Belanger
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

#include <PapillonNDL/angle_table.hpp>
#include <PapillonNDL/equiprobable_angle_bins.hpp>
#include <PapillonNDL/isotropic.hpp>

namespace py = pybind11;

using namespace pndl;

// Trampoline class for abstract pndl::AngleLaw
class PyAngleLaw : public AngleLaw {
 public:
  using AngleLaw::AngleLaw;

  double sample_mu(double xi) const override {
    PYBIND11_OVERRIDE_PURE(double, AngleLaw, sample_mu, xi);
  }
};

void init_AngleLaw(py::module& m) {
  py::class_<AngleLaw, PyAngleLaw, std::shared_ptr<AngleLaw>>(m, "AngleLaw")
      .def(py::init<>())
      .def("sample_mu", &AngleLaw::sample_mu);
}

void init_Isotropic(py::module& m) {
  py::class_<Isotropic, AngleLaw, std::shared_ptr<Isotropic>>(m, "Isotropic")
      .def(py::init<>())
      .def("sample_mu", &Isotropic::sample_mu);
}

void init_EquiprobableAngleBins(py::module& m) {
  py::class_<EquiprobableAngleBins, AngleLaw,
             std::shared_ptr<EquiprobableAngleBins>>(m, "EquiprobableAngleBins")
      .def(py::init<const ACE&, size_t>())
      .def("sample_mu", &EquiprobableAngleBins::sample_mu)
      .def("size", &EquiprobableAngleBins::size)
      .def("bin_bounds", &EquiprobableAngleBins::bin_bounds);
}

void init_AngleTable(py::module& m) {
  py::class_<AngleTable, AngleLaw, std::shared_ptr<AngleTable>>(m, "AngleTable")
      .def(py::init<const ACE&, size_t>())
      .def("sample_mu", &AngleTable::sample_mu)
      .def("size", &AngleTable::size)
      .def("cosines", &AngleTable::cosines)
      .def("pdf", &AngleTable::pdf)
      .def("cdf", &AngleTable::cdf)
      .def("interpolate", &AngleTable::interpolation);
}
