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

#include <PapillonNDL/multi_region_1d.hpp>

namespace py = pybind11;

using namespace pndl;

void init_MultiRegion1D(py::module& m) {
  py::class_<MultiRegion1D, Tabulated1D, std::shared_ptr<MultiRegion1D>>(m, "MultiRegion1D")
    .def(py::init<const std::vector<Region1D>&>())
    .def(py::init<const std::vector<uint32_t>&, const std::vector<Interpolation>&, const std::vector<double>&, const std::vector<double>&>())
    .def("__call__", &MultiRegion1D::operator())
    .def("integrate", &MultiRegion1D::integrate)
    .def("breakpoints", &MultiRegion1D::breakpoints)
    .def("interpolation", &MultiRegion1D::interpolation)
    .def("x", &MultiRegion1D::x)
    .def("y", &MultiRegion1D::y)
    .def("size", &MultiRegion1D::size)
    .def("min_x", &MultiRegion1D::min_x)
    .def("max_x", &MultiRegion1D::max_x)
  ;
}
