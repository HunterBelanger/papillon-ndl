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

#include <PapillonNDL/energy_grid.hpp>
#include <memory>

namespace py = pybind11;

using namespace pndl;

void init_EnergyGrid(py::module& m) {
  py::class_<EnergyGrid, std::shared_ptr<EnergyGrid>>(m, "EnergyGrid")
      .def(py::init<const ACE&, uint32_t>())
      .def(py::init<const std::vector<double>&, uint32_t>())
      .def("__getitem__", &EnergyGrid::operator[])
      .def("grid", &EnergyGrid::grid)
      .def("size", &EnergyGrid::size)
      .def("min_energy", &EnergyGrid::min_energy)
      .def("max_energy", &EnergyGrid::max_energy)
      .def("get_lower_index", &EnergyGrid::get_lower_index)
      .def("hash_energy_grid", &EnergyGrid::hash_energy_grid);
}
