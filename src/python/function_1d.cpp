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

#include <PapillonNDL/constant.hpp>
#include <PapillonNDL/multi_region_1d.hpp>
#include <PapillonNDL/polynomial_1d.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <PapillonNDL/sum_1d.hpp>
#include <PapillonNDL/difference_1d.hpp>

namespace py = pybind11;

using namespace pndl;

// Trampoline class for abstract pndl::Function1D
class PyFunction1D : public Function1D {
 public:
  using Function1D::Function1D;

  double operator()(double x) const override {
    PYBIND11_OVERRIDE_PURE_NAME(double, Function1D, "__call__", operator(), x);
  }

  double integrate(double x_low, double x_hi) const override {
    PYBIND11_OVERRIDE_PURE(double, Function1D, integrate, x_low, x_hi);
  }
};

// Trampoline class for abstract pndl::Tabulated1D
class PyTabulated1D : public Tabulated1D {
 public:
  using Tabulated1D::Tabulated1D;

  double operator()(double x) const override {
    PYBIND11_OVERRIDE_PURE_NAME(double, Tabulated1D, "__call__", operator(), x);
  }

  double integrate(double x_low, double x_hi) const override {
    PYBIND11_OVERRIDE_PURE(double, Tabulated1D, integrate, x_low, x_hi);
  }

  std::vector<uint32_t> breakpoints() const override {
    PYBIND11_OVERRIDE_PURE(std::vector<uint32_t>, Tabulated1D, breakpoints);
  }

  std::vector<Interpolation> interpolation() const override {
    PYBIND11_OVERRIDE_PURE(std::vector<Interpolation>, Tabulated1D,
                           interpolation);
  }

  std::vector<double> x() const override {
    PYBIND11_OVERRIDE_PURE(std::vector<double>, Tabulated1D, x);
  }

  std::vector<double> y() const override {
    PYBIND11_OVERRIDE_PURE(std::vector<double>, Tabulated1D, y);
  }
};

void init_Function1D(py::module& m) {
  py::class_<Function1D, PyFunction1D, std::shared_ptr<Function1D>>(
      m, "Function1D")
      .def(py::init<>())
      .def("__call__", &Function1D::operator())
      .def("evaluate", &Function1D::evaluate)
      .def("integrate", &Function1D::integrate);
}

void init_Constant(py::module& m) {
  py::class_<Constant, Function1D, std::shared_ptr<Constant>>(m, "Constant")
      .def(py::init<double>())
      .def("__call__", &Constant::operator())
      .def("evaluate", &Constant::evaluate)
      .def("integrate", &Constant::integrate);
}

void init_Polynomial1D(py::module& m) {
  py::class_<Polynomial1D, Function1D, std::shared_ptr<Polynomial1D>>(
      m, "Polynomial1D")
      .def(py::init<const std::vector<double>&>())
      .def("__call__", &Polynomial1D::operator())
      .def("evaluate", &Polynomial1D::evaluate)
      .def("integrate", &Polynomial1D::integrate)
      .def("order", &Polynomial1D::order)
      .def("coefficient", &Polynomial1D::coefficient);
}

void init_Tabulated1D(py::module& m) {
  py::class_<Tabulated1D, PyTabulated1D, std::shared_ptr<Tabulated1D>>(
      m, "Tabulated1D")
      .def(py::init<>())
      .def("__call__", &Tabulated1D::operator())
      .def("evaluate", &Tabulated1D::evaluate)
      .def("breakpoints", &Tabulated1D::breakpoints)
      .def("interpolation", &Tabulated1D::interpolation)
      .def("x", &Tabulated1D::x)
      .def("y", &Tabulated1D::y);
}

void init_Region1D(py::module& m) {
  py::class_<Region1D, Tabulated1D, std::shared_ptr<Region1D>>(m, "Region1D")
      .def(py::init<const std::vector<double>&, const std::vector<double>&,
                    Interpolation>())
      .def("__call__", &Region1D::operator())
      .def("evaluate", &Region1D::evaluate)
      .def("integrate", &Region1D::integrate)
      .def("breakpoints", &Region1D::breakpoints)
      .def("interpolation", &Region1D::interpolation)
      .def("x", &Region1D::x)
      .def("y", &Region1D::y)
      .def("size", &Region1D::size)
      .def("min_x", &Region1D::min_x)
      .def("max_x", &Region1D::max_x);
}

void init_MultiRegion1D(py::module& m) {
  py::class_<MultiRegion1D, Tabulated1D, std::shared_ptr<MultiRegion1D>>(
      m, "MultiRegion1D")
      .def(py::init<const std::vector<Region1D>&>())
      .def(py::init<const std::vector<uint32_t>&,
                    const std::vector<Interpolation>&,
                    const std::vector<double>&, const std::vector<double>&>())
      .def("__call__", &MultiRegion1D::operator())
      .def("evaluate", &MultiRegion1D::evaluate)
      .def("integrate", &MultiRegion1D::integrate)
      .def("breakpoints", &MultiRegion1D::breakpoints)
      .def("interpolation", &MultiRegion1D::interpolation)
      .def("x", &MultiRegion1D::x)
      .def("y", &MultiRegion1D::y)
      .def("size", &MultiRegion1D::size)
      .def("min_x", &MultiRegion1D::min_x)
      .def("max_x", &MultiRegion1D::max_x);
}

void init_Sum1D(py::module& m) {
  py::class_<Sum1D, Function1D, std::shared_ptr<Sum1D>>(
      m, "Sum1D")
      .def(py::init<std::shared_ptr<Function1D>, std::shared_ptr<Function1D>>())
      .def("__call__", &Sum1D::operator())
      .def("evaluate", &Sum1D::evaluate)
      .def("integrate", &Sum1D::integrate)
      .def("term_1", &Sum1D::term_1)
      .def("term_2", &Sum1D::term_2);
}

void init_Difference1D(py::module& m) {
  py::class_<Difference1D, Function1D, std::shared_ptr<Difference1D>>(
      m, "Difference1D")
      .def(py::init<std::shared_ptr<Function1D>, std::shared_ptr<Function1D>>())
      .def("__call__", &Difference1D::operator())
      .def("evaluate", &Difference1D::evaluate)
      .def("integrate", &Difference1D::integrate)
      .def("term_1", &Difference1D::term_1)
      .def("term_2", &Difference1D::term_2);
}