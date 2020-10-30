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
#include <pybind11/functional.h>

namespace py = pybind11;

// LCG parameters
constexpr uint64_t g   {2806196910506780709LL};   // multiplication
constexpr uint64_t c    {1};                       // additive factor, c
constexpr uint64_t M    {0x8000000000000000};      // 2^63
constexpr uint64_t mask   {0x7fffffffffffffff};      // 2^63 - 1
constexpr double   normalization   {1.0 / M};           // 2^-63
uint64_t rng_seed_ = 1;

void seed(uint64_t s) {
  rng_seed_ = s;
}

uint32_t get_current_seed() {
  return rng_seed_;
}

void advance_seed(uint64_t n) {
  n &= mask;

  uint64_t g_tmp = g;
  uint64_t c_tmp = c;
  uint64_t g_new = 1;
  uint64_t c_new = 0;

  while (n > 0) {
    if (n & 1) {
      g_new *= g_tmp;
      c_new = c_new * g_tmp + c_tmp;
    }
    c_tmp *= (g_tmp + 1);
    g_tmp *= g_tmp;

    n >>= 1;
  }

  rng_seed_ = (g_new * rng_seed_ + c_new) & mask;
}

double rang() {
  rng_seed_ = (g*rng_seed_ + c) & mask;

  return rng_seed_ * normalization;
}

std::function<double()> rng(rang);

void init_PRNG(py::module& m) {
  m.def("seed", &seed);
  m.def("get_current_seed", &get_current_seed);
  m.def("advance_seed", &advance_seed);
  m.def("rang", &rang);
  m.attr("rng") = rng;
}
