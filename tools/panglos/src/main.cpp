/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2022, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Papillon Nuclear Data Library (PapillonNDL).
 *
 * PapillonNDL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PapillonNDL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PapillonNDL. If not, see <https://www.gnu.org/licenses/>.
 *
 * */

/**
 * @file
 * @author Hunter Belanger
 */

#include <cstddef>
#include <iostream>
#include <string>
#include <memory>
#include <variant>

#include <ndarray.hpp>
#include <boost/hana.hpp>  // Needed for the _c literal for constructing mt4
#include <ENDFtk.hpp>
using namespace njoy::ENDFtk;

#include "incoherent_inelastic.hpp"
#include "incoherent_elastic.hpp"
#include "coherent_elastic.hpp"

int main(const int argc, const char** argv) {
  //std::string fname = "tsl-HinH2O.endf";
  //int MAT = 1;

  std::string fname = "tsl-reactor-graphite-10P.endf";
  int MAT = 31;

  tree::Tape<std::string> pendf = tree::fromFile(fname);
  file::Type<7> mf7 = pendf.material(MAT).front().file(7).parse<7>();

  section::Type<7, 4> mt4 = mf7.section(4_c); 
  IncoherentInelastic ii(mt4);
  const double Emin = ii.Emin();
  const double Emax = ii.Emax();
  const std::size_t ti = 0;

  // Make energy grid in log scale
  {
    constexpr std::size_t NE = 100;
    std::vector<double> Egrid(NE, 0.);
    const double du =
        (std::log(Emax) - std::log(Emin)) / (static_cast<double>(NE) - 1.);
    for (std::size_t i = 0; i < NE; i++) {
      Egrid[i] = std::exp(std::log(Emin) + static_cast<double>(i) * du);
    }
    if (Egrid.back() > Emax) Egrid.back() = Emax;

    NDArray<double> IIxs({2, NE});
    for (std::size_t i = 0; i < NE; i++) {
      const auto& E = Egrid[i];
      std::cout << "Running index " << i << ", Energy " << E << " eV.\n";
      IIxs(0, i) = E;
      IIxs(1, i) = ii.xs(ti, E);
    }
    IIxs.save("iixs.npy");
  }

  // Do Elastic
  std::unique_ptr<CoherentElastic> ce = nullptr;
  std::unique_ptr<IncoherentElastic> ie = nullptr;
  if (mf7.hasSection(2)) {
    section::Type<7, 2> mt2 = mf7.section(2_c);  
    const auto& scatter_law = mt2.scatteringLaw();

    if (std::holds_alternative<section::Type<7,2>::CoherentElastic>(scatter_law)) {
      ce = std::make_unique<CoherentElastic>(std::get<section::Type<7,2>::CoherentElastic>(scatter_law));
    } else if (std::holds_alternative<section::Type<7,2>::IncoherentElastic>(scatter_law)) {
      ie = std::make_unique<IncoherentElastic>(std::get<section::Type<7,2>::IncoherentElastic>(scatter_law));
    } else {
      const section::Type<7,2>::MixedElastic& me = std::get<section::Type<7,2>::MixedElastic>(scatter_law);
      ce = std::make_unique<CoherentElastic>(me.coherent());
      ie = std::make_unique<IncoherentElastic>(me.incoherent());
    }
  }

  if (ce) {
    constexpr std::size_t NE = 5000;
    std::vector<double> Egrid(NE, 0.);
    const double du =
        (std::log(Emax) - std::log(Emin)) / (static_cast<double>(NE) - 1.);
    for (std::size_t i = 0; i < NE; i++) {
      Egrid[i] = std::exp(std::log(Emin) + static_cast<double>(i) * du);
    }
    if (Egrid.back() > Emax) Egrid.back() = Emax;

    const double T = ii.temperatures()[ti];
    NDArray<double> CExs({2, NE});
    for (std::size_t i = 0; i < NE; i++) {
      const auto& E = Egrid[i];
      CExs(0, i) = E;
      CExs(1, i) = ce->xs(T, E);
    }
    CExs.save("cexs.npy");
  }
  
  return 0;
}
